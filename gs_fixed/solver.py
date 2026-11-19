"""Fixed-boundary P3 GS, Picard, and prescribed-q solvers.

Parallelism is between independent F/FFprime candidate equilibria.
Each process owns a persistent sparse factorization and field locator.
"""
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from dataclasses import dataclass, field, fields, replace
import multiprocessing as mp
import os
import time
import warnings
from functools import wraps
from threadpoolctl import threadpool_limits
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree
from .mesh import P1Mesh, P3Mesh, build_p1_mesh, generate_p3_nodes
from .fem import (precompute_quadrature, assemble_k_matrix, assemble_rhs,
                  apply_bc, _prepare_rhs_profiles)
from .profiles import (FFprimeModel, evaluate_ffprime_model,
    initialize_ffprime_model, normalize_q_target, evaluate_q_target,
    reconstruct_g_from_ffprime, reconstruct_pressure, validate_profile_grid)
from .post import (QOptions,P3FieldEvaluator,build_q_contour_cache,
                    evaluate_q_geometry,evaluate_q_from_g)


def _single_blas_thread(function):
    @wraps(function)
    def controlled(*args,**kwargs):
        with threadpool_limits(limits=1):
            return function(*args,**kwargs)
    return controlled


def _options(cls,value):
    if value is None: return cls()
    if isinstance(value,cls): return value
    if isinstance(value,dict): return cls(**value)
    raise TypeError(f"Expected {cls.__name__} or a dictionary.")


def _vector(value,name,size=None):
    raw=np.asarray(value)
    if np.iscomplexobj(raw): raise ValueError(f"{name} must be real.")
    a=np.array(raw,dtype=float,copy=True)
    if a.ndim not in (1,2) or (a.ndim==2 and 1 not in a.shape):
        raise ValueError(f"{name} must be a vector.")
    a=a.reshape(-1)
    if not np.all(np.isfinite(a)) or (size is not None and len(a)!=size):
        raise ValueError(f"Invalid {name}.")
    return a


def _positive(value,name,allow_zero=False):
    if isinstance(value,(bool,np.bool_)) or np.ndim(value)!=0 or np.iscomplexobj(value):
        raise ValueError(f"Invalid {name}.")
    x=float(value)
    if not np.isfinite(x) or x<0 or (x==0 and not allow_zero):
        raise ValueError(f"Invalid {name}.")
    return x


def _integer(value,name,lower=1):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or value<lower:
        raise ValueError(f"{name} must be an integer >= {lower}.")
    return int(value)


@dataclass(frozen=True)
class PicardOptions:
    omega: float = .5
    max_iterations: int = 100
    update_tolerance: float = 1e-8
    residual_tolerance: float = 1e-8
    axis_mode: str = "max"
    symmetry_tolerance: float = 1e-12
    flux_span_tolerance: float = 1e-12
    verbose: bool = False
    error_on_failure: bool = False

    def __post_init__(self):
        if _positive(self.omega,"omega")>1: raise ValueError("omega must be <=1.")
        _integer(self.max_iterations,"max_iterations")
        for name in ("update_tolerance","residual_tolerance","symmetry_tolerance","flux_span_tolerance"):
            _positive(getattr(self,name),name)
        if self.axis_mode not in ("max","min","auto"): raise ValueError("Invalid axis_mode.")


@dataclass(frozen=True)
class PQOptions:
    max_iterations: int = 20
    minimum_iterations: int = 1
    q_relative_tolerance: float = 1e-3
    q_rms_log_tolerance: float = 1e-3
    minimum_q_points: int = 3
    fit_range: tuple = (0.,.95)
    axis_weight: float = .25
    initial_scales: tuple = (1.,.5,0.)
    finite_difference_absolute_step: float = .05
    finite_difference_relative_step: float = .02
    regularization: float = 1e-4
    jacobian_ridge: float = 1e-10
    initial_damping: float = 1e-2
    damping_growth: float = 10.
    maximum_damping_attempts: int = 4
    maximum_line_search_steps: int = 6
    max_coefficient_step: float = 2.
    minimum_coefficient_step: float = 1e-8
    armijo_fraction: float = 1e-3
    minimum_absolute_merit_decrease: float = 1e-12
    minimum_relative_merit_decrease: float = 1e-8
    error_on_failure: bool = False
    verbose: bool = True

    def __post_init__(self):
        for name in ("max_iterations","minimum_iterations","minimum_q_points",
                     "maximum_damping_attempts","maximum_line_search_steps"):
            _integer(getattr(self,name),name)
        if self.minimum_iterations>self.max_iterations: raise ValueError("minimum_iterations exceeds maximum.")
        for name in ("q_relative_tolerance","q_rms_log_tolerance","axis_weight",
            "finite_difference_absolute_step","finite_difference_relative_step",
            "initial_damping","max_coefficient_step","minimum_coefficient_step"):
            _positive(getattr(self,name),name)
        for name in ("regularization","jacobian_ridge","armijo_fraction",
                     "minimum_absolute_merit_decrease","minimum_relative_merit_decrease"):
            _positive(getattr(self,name),name,True)
        if _positive(self.damping_growth,"damping_growth")<=1: raise ValueError("damping_growth must exceed 1.")
        f=np.asarray(self.fit_range,dtype=float)
        if f.shape!=(2,) or not np.all(np.isfinite(f)) or not 0<=f[0]<f[1]<=1:
            raise ValueError("fit_range must be an increasing interval in [0,1].")
        scales=np.asarray(self.initial_scales,dtype=float)
        if scales.ndim!=1 or not len(scales) or not np.all(np.isfinite(scales)) or np.any(scales<0):
            raise ValueError("initial_scales must be finite and nonnegative.")
        object.__setattr__(self,"fit_range",tuple(f))
        object.__setattr__(self,"initial_scales",tuple(scales))


@dataclass(frozen=True)
class ParallelOptions:
    max_workers: int = 1
    blas_threads: int = 1
    start_method: str = "spawn"

    def __post_init__(self):
        _integer(self.max_workers,"max_workers",0) # 0: allocation-aware automatic.
        _integer(self.blas_threads,"blas_threads")
        if self.start_method not in mp.get_all_start_methods():
            raise ValueError("Unsupported multiprocessing start method.")


@dataclass(frozen=True)
class RHSCache:
    operator: object
    radius: np.ndarray
    weight: np.ndarray

    def assemble(self,psi,pprime,ffprime,dpsi):
        span,profiles=_prepare_rhs_profiles(pprime,ffprime,dpsi)
        psi_n=(np.asarray(self.operator@psi)+span)/span
        if not np.all(np.isfinite(psi_n)): raise ValueError("Nonfinite normalized flux.")
        physical=profiles(psi_n)
        source=4*np.pi*1e-7*self.radius*physical[:,0]+physical[:,1]/self.radius
        rhs=np.asarray(self.operator.T@(self.weight*source)).reshape(-1)
        if not np.all(np.isfinite(rhs)): raise ValueError("Nonfinite RHS.")
        return rhs


def build_rhs_cache(mesh,quad,max_mb=512.):
    """Sparse local-to-quadrature operator: same integration as assemble_rhs."""
    max_mb=_positive(max_mb,"max_mb")
    nt=len(mesh.elements); nq=quad.nq
    entries=nt*nq*10
    # Bound triplet construction plus the final CSR and coordinate arrays.
    estimated=entries*40+(nt*nq+1)*16
    if estimated>max_mb*1024**2: return None
    rows=np.repeat(np.arange(nt*nq,dtype=np.int64),10)
    columns=np.repeat(mesh.elements,nq,axis=0).ravel()
    values=np.tile(quad.N.ravel(),nt)
    op=sparse.coo_matrix((values,(rows,columns)),shape=(nt*nq,mesh.n_dof)).tocsr()
    vertices=mesh.p1_points[mesh.p1_elements]
    e1=vertices[:,1]-vertices[:,0]; e2=vertices[:,2]-vertices[:,0]
    det=e1[:,0]*e2[:,1]-e2[:,0]*e1[:,1]
    radius=(vertices[:,0,0,None]+e1[:,0,None]*quad.xi+e2[:,0,None]*quad.eta).ravel()
    weight=(det[:,None]*quad.w).ravel()
    if np.any(radius<=0): raise ValueError("Nonpositive quadrature radius.")
    return RHSCache(op,radius,weight)


@dataclass
class PreparedGSProblem:
    mesh: P3Mesh
    psi_n: np.ndarray
    r0: float
    b0: float
    a_minor: float
    picard: PicardOptions
    q_options: QOptions
    quad: object
    k: object
    kii: object
    bc: object
    initial_psi: np.ndarray
    rhs_cache: object
    symmetry_error: float
    preparation_seconds: float
    factor: object = field(repr=False)

    def assemble_rhs(self,psi,pprime,ffprime,dpsi):
        if self.rhs_cache is not None:
            return self.rhs_cache.assemble(psi,pprime,ffprime,dpsi)
        return assemble_rhs(self.mesh,psi,self.quad,pprime,ffprime,dpsi)

    def worker_payload(self):
        # SuperLU itself is neither sent nor shared between processes.
        return {f.name:getattr(self,f.name) for f in fields(self) if f.name!="factor"}

    @classmethod
    def from_payload(cls,payload):
        return cls(**payload,factor=splu(payload["kii"].tocsc()))


@_single_blas_thread
def prepare_gs_problem(boundary_points=None,*,mesh=None,r0,b0,a_minor,
                       nr=55,nz=65,n_profile=201,quadrature_degree=10,
                       picard=None,q_options=None,psi0=None,psi_axis_guess=None,
                       cache_rhs=True,rhs_cache_max_mb=512.):
    """Prepare and factor a zero-boundary-flux problem once per geometry."""
    start=time.perf_counter()
    r0=_positive(r0,"r0"); a_minor=_positive(a_minor,"a_minor")
    if np.ndim(b0)!=0 or np.iscomplexobj(b0) or not np.isfinite(b0) or float(b0)==0:
        raise ValueError("b0 must be a finite nonzero scalar.")
    b0=float(b0)
    _integer(n_profile,"n_profile",2)
    options=_options(PicardOptions,picard); q_options=_options(QOptions,q_options)
    if mesh is None:
        if boundary_points is None: raise ValueError("Provide boundary_points or mesh.")
        p1=build_p1_mesh(boundary_points,nr=nr,nz=nz,extra_points=[[r0,0.]])
        mesh=generate_p3_nodes(p1)
    elif isinstance(mesh,P1Mesh): mesh=generate_p3_nodes(mesh)
    if not isinstance(mesh,P3Mesh): raise TypeError("mesh must be P1Mesh or P3Mesh.")
    s=validate_profile_grid(np.linspace(0,1,n_profile),require_uniform=True)
    quad=precompute_quadrature(quadrature_degree)
    k=assemble_k_matrix(mesh,quad)
    kii,_,bc=apply_bc(mesh,k,np.zeros(mesh.n_dof),0.)
    symmetry=sparse.linalg.norm(kii-kii.T)/max(sparse.linalg.norm(kii),np.finfo(float).tiny)
    if symmetry>options.symmetry_tolerance: raise ValueError("Stiffness matrix is not symmetric.")
    factor=splu(kii.tocsc())
    if psi0 is None:
        guess=(-.5 if options.axis_mode=="min" else .5) if psi_axis_guess is None else float(psi_axis_guess)
        if not np.isfinite(guess) or abs(guess)<=100*np.spacing(max(1.,abs(guess))):
            raise ValueError("Initial axis flux must differ from zero.")
        boundary=np.asarray(boundary_points) if boundary_points is not None else mesh.p1_points[np.unique(mesh.boundary_edges)]
        distance=cKDTree(boundary).query(mesh.points)[0]
        if distance.max()<=0: raise ValueError("Cannot construct initial flux.")
        initial=guess*distance/distance.max()
    else: initial=_vector(psi0,"psi0",mesh.n_dof)
    initial[bc.boundary_nodes]=0.
    rhs_cache=build_rhs_cache(mesh,quad,rhs_cache_max_mb) if cache_rhs else None
    return PreparedGSProblem(mesh,s,r0,b0,a_minor,options,q_options,quad,k,kii,bc,
                             initial,rhs_cache,float(symmetry),time.perf_counter()-start,factor)


@dataclass(frozen=True)
class PicardResult:
    converged: bool
    iterations: int
    update_history: np.ndarray
    residual_history: np.ndarray
    psi_axis_history: np.ndarray
    dpsi_history: np.ndarray
    min_psi_n_history: np.ndarray
    max_psi_n_history: np.ndarray
    initial_boundary_error: float
    elapsed_seconds: float


@dataclass(frozen=True)
class GSState:
    psi: np.ndarray
    psi_n: np.ndarray
    psi_axis: float
    dpsi: float
    axis_mode: str
    axis_node: int
    axis_point: np.ndarray
    final_rhs: np.ndarray
    final_residual: np.ndarray
    inner_converged: bool
    inner_iterations: int


def _axis_and_span(psi,free,mode,tolerance):
    values=psi[free]
    if mode=="auto":
        mode="max" if abs(values.max())>=abs(values.min()) else "min"
    local=int(np.argmax(values) if mode=="max" else np.argmin(values))
    node=int(free[local]); axis=float(psi[node]); span=-axis
    scale=max(1.,abs(axis),np.linalg.norm(psi,np.inf))
    if abs(span)<=tolerance*scale: raise ValueError("Axis and boundary flux are indistinguishable.")
    return axis,span,mode,node


def solve_gs_inner(problem,pprime,ffprime,psi0=None,options=None,*,emit_warning=True):
    """Picard iteration with a reused factorization and true nonlinear residual."""
    opts=problem.picard if options is None else _options(PicardOptions,options)
    if problem.symmetry_error>opts.symmetry_tolerance: raise ValueError("Stiffness symmetry tolerance exceeded.")
    pp=_vector(pprime,"pprime",len(problem.psi_n))
    ff=_vector(ffprime,"ffprime",len(problem.psi_n))
    psi=_vector(problem.initial_psi if psi0 is None else psi0,"psi0",problem.mesh.n_dof)
    start=time.perf_counter(); free=problem.bc.free_nodes; boundary=problem.bc.boundary_nodes
    initial_error=float(np.max(abs(psi[boundary]))); psi[boundary]=0.
    axis,span,mode,node=_axis_and_span(psi,free,opts.axis_mode,opts.flux_span_tolerance)
    rhs=problem.assemble_rhs(psi,pp,ff,span)
    histories=[]; converged=False
    for iteration in range(1,opts.max_iterations+1):
        solved=problem.factor.solve(problem.bc.reduce_rhs(rhs))
        new=psi.copy()
        new[free]=(1-opts.omega)*psi[free]+opts.omega*solved
        new[boundary]=0.
        update=np.linalg.norm(new[free]-psi[free])/max(np.linalg.norm(new[free]),
                                                      np.linalg.norm(psi[free]),np.finfo(float).eps)
        axis,span,mode,node=_axis_and_span(new,free,opts.axis_mode,opts.flux_span_tolerance)
        next_rhs=problem.assemble_rhs(new,pp,ff,span)
        lhs=np.asarray(problem.k[free,:]@new)
        residual=lhs-next_rhs[free]
        relative=np.linalg.norm(residual)/max(np.linalg.norm(lhs),np.linalg.norm(next_rhs[free]),np.finfo(float).eps)
        psi_n=(new-axis)/span
        histories.append((update,relative,axis,span,psi_n.min(),psi_n.max()))
        psi,rhs=new,next_rhs
        if opts.verbose: print(f"Picard {iteration:3d}: update={update:.3e}, residual={relative:.3e}",flush=True)
        if update<opts.update_tolerance and relative<opts.residual_tolerance:
            converged=True; break
    h=np.asarray(histories)
    result=PicardResult(converged,iteration,*[h[:,j] for j in range(6)],
                         initial_error,time.perf_counter()-start)
    state=GSState(psi,(psi-axis)/span,axis,span,mode,node,problem.mesh.points[node].copy(),
                   rhs,residual,converged,iteration)
    if not converged:
        message=f"GS Picard failed after {iteration} iterations: update={h[-1,0]:.3e}, residual={h[-1,1]:.3e}."
        if opts.error_on_failure: raise RuntimeError(message)
        if emit_warning: warnings.warn(message,RuntimeWarning,stacklevel=2)
    return state,result



@dataclass(frozen=True)
class QMetrics:
    valid_mask: np.ndarray
    relative_residual: np.ndarray
    log_residual: np.ndarray
    weights: np.ndarray
    max_abs_relative: float
    rms_log: float
    weighted_merit: float


def calculate_q_metrics(q_result,target,options):
    s,q,factor=q_result.psi_n,q_result.q,q_result.q_per_f
    valid=(np.isfinite(s)&np.isfinite(q)&np.isfinite(target)&np.isfinite(factor)
           &(q!=0)&(target!=0)&(factor!=0)&(s>=options.fit_range[0])&(s<=options.fit_range[1]))
    if np.any(valid&(np.sign(q)!=np.sign(target))): raise ValueError("Calculated and prescribed q signs disagree.")
    if np.count_nonzero(valid)<options.minimum_q_points: raise ValueError("Too few valid q surfaces.")
    relative=np.full_like(q,np.nan); log=np.full_like(q,np.nan); weights=np.zeros_like(q)
    relative[valid]=(q[valid]-target[valid])/target[valid]
    log[valid]=np.log(np.abs(q[valid]/target[valid]))
    weights[valid]=1.
    weights[valid&(np.abs(s)<=100*np.finfo(float).eps)]=options.axis_weight
    merit=.5*np.sum(weights[valid]*log[valid]**2)/weights[valid].sum()
    return QMetrics(valid,relative,log,weights,float(np.max(abs(relative[valid]))),
                     float(np.sqrt(np.mean(log[valid]**2))),float(merit))


@dataclass
class _Candidate:
    coefficients: np.ndarray
    valid: bool = False
    failure: str = ""
    state: object = None
    inner: object = None
    profile: object = None
    q_result: object = None
    metrics: object = None
    pid: int = field(default_factory=os.getpid)

    @property
    def residual(self): return self.metrics.log_residual[self.metrics.valid_mask]
    @property
    def weights(self): return self.metrics.weights[self.metrics.valid_mask]
    @property
    def merit(self): return self.metrics.weighted_merit if self.valid else np.inf


class _Runtime:
    def __init__(self,problem,pprime,model,target,inner,outer,cache=None):
        self.problem=problem; self.pprime=pprime; self.model=model; self.target=target
        self.inner=replace(inner,verbose=False,error_on_failure=False); self.outer=outer
        self.evaluator=P3FieldEvaluator(problem.mesh)
        self.cache=(build_q_contour_cache(problem.mesh,problem.q_options,self.evaluator)
                    if cache is None else cache)

    def candidate(self,job):
        coefficients,psi_guess,required_mode,required_sign=job
        candidate=_Candidate(np.asarray(coefficients).copy())
        try:
            profile=evaluate_ffprime_model(self.model.with_coefficients(coefficients),self.problem.psi_n)
            state,inner=solve_gs_inner(self.problem,self.pprime,profile.ffprime,psi_guess,
                                      self.inner,emit_warning=False)
            candidate.state=state; candidate.inner=inner
            if not inner.converged: raise ValueError(f"inner_GS_nonconvergence: residual={inner.residual_history[-1]:.3e}")
            if required_mode and state.axis_mode!=required_mode: raise ValueError("magnetic_axis_branch_change")
            if required_sign and np.sign(state.dpsi)!=required_sign: raise ValueError("flux_span_sign_change")
            geometry=evaluate_q_geometry(self.problem.mesh,state,self.problem.q_options,
                                         evaluator=self.evaluator,cache=self.cache)
            q_result=evaluate_q_from_g(geometry,profile)
            metrics=calculate_q_metrics(q_result,self.target,self.outer)
            candidate.profile=profile; candidate.q_result=q_result
            candidate.metrics=metrics; candidate.valid=True
        except Exception as exc:
            candidate.failure=f"{type(exc).__name__}: {exc}"
        return candidate

    def difference(self,job):
        j,coefficients,psi,residual,mask,mode,sign=job
        step=max(self.outer.finite_difference_absolute_step,
                 self.outer.finite_difference_relative_step*max(1.,abs(coefficients[j])))
        failures=[]
        for n,direction in enumerate((1.,-1.),start=1):
            trial=coefficients.copy(); trial[j]+=direction*step
            candidate=self.candidate((trial,psi,mode,sign))
            if candidate.valid and np.array_equal(candidate.metrics.valid_mask,mask):
                column=(candidate.residual-residual)/(direction*step)
                return column,dict(success=True,step=step,direction=direction,evaluations=n,pid=os.getpid())
            failures.append(candidate.failure or "changed q residual mask")
        return None,dict(success=False,step=step,direction=0.,evaluations=2,pid=os.getpid(),failures=failures)


_WORKER_RUNTIME=None
_WORKER_THREAD_LIMIT=None


def _worker_init(payload,pprime,model,target,inner,outer,threads,cache):
    global _WORKER_RUNTIME,_WORKER_THREAD_LIMIT
    from threadpoolctl import threadpool_limits
    _WORKER_THREAD_LIMIT=threadpool_limits(limits=threads)
    problem=PreparedGSProblem.from_payload(payload)
    _WORKER_RUNTIME=_Runtime(problem,pprime,model,target,inner,outer,cache)


def _worker_candidate(job): return _WORKER_RUNTIME.candidate(job)
def _worker_difference(job): return _WORKER_RUNTIME.difference(job)


def available_cpus():
    limits=[os.cpu_count() or 1]
    if hasattr(os,"sched_getaffinity"):
        limits.append(len(os.sched_getaffinity(0)))
    value=os.environ.get("SLURM_CPUS_PER_TASK")
    if value:
        try:
            if int(value)>0: limits.append(int(value))
        except ValueError: pass
    return max(1,min(limits))


class _CandidatePool:
    """Lazily start a persistent process pool only for concurrent work."""
    def __init__(self,runtime,parallel,max_jobs):
        self.runtime=runtime; self.parallel=parallel; self.executor=None
        allocated=max(1,available_cpus()//parallel.blas_threads)
        requested=parallel.max_workers or allocated
        self.workers=max(1,min(requested,allocated,max_jobs))
        self.pids=set()

    def map(self,kind,jobs):
        if self.workers<=1 or len(jobs)<=1:
            result=[getattr(self.runtime,kind)(job) for job in jobs]
        else:
            if self.executor is None:
                r=self.runtime
                self.executor=ProcessPoolExecutor(
                    max_workers=self.workers,mp_context=mp.get_context(self.parallel.start_method),
                    initializer=_worker_init,initargs=(r.problem.worker_payload(),r.pprime,
                        r.model,r.target,r.inner,r.outer,self.parallel.blas_threads,r.cache))
            function=_worker_candidate if kind=="candidate" else _worker_difference
            result=list(self.executor.map(function,jobs,chunksize=1))
        for value in result:
            self.pids.add(value.pid if kind=="candidate" else value[1]["pid"])
        return result

    def close(self):
        if self.executor is not None: self.executor.shutdown(wait=True)


def _q_converged(candidate,options):
    return (candidate.metrics.max_abs_relative<=options.q_relative_tolerance
            and candidate.metrics.rms_log<=options.q_rms_log_tolerance)


def _lm_step(candidate,response,damping,options):
    root=np.sqrt(candidate.weights)
    jacobian=response*root[:,None]; residual=candidate.residual*root
    h=jacobian.T@jacobian
    n=len(candidate.coefficients)
    regularizer=np.diag(np.arange(n,dtype=float)**2)
    scale=max(np.trace(h)/n,1e-10)
    normal=(h+options.regularization*regularizer
            +damping*np.diag(np.maximum(np.diag(h),scale))
            +options.jacobian_ridge*np.eye(n))
    gradient=jacobian.T@residual+options.regularization*regularizer@candidate.coefficients
    step=-np.linalg.solve(normal,gradient)
    predicted=residual+jacobian@step
    reduction=.5*(residual@residual-predicted@predicted)/max(candidate.weights.sum(),np.finfo(float).eps)
    return step,float(reduction)


@dataclass(frozen=True)
class PQSolveResult:
    converged: bool
    termination_reason: str
    outer_iterations: int
    function_evaluations: int
    elapsed_seconds: float
    final_model: object
    final_f_profile: object
    final_inner_result: PicardResult
    final_q_result: object
    final_metrics: QMetrics
    history: tuple
    initial_trials: tuple
    last_update: dict
    parallel_info: dict


def solve_pq_constraint(problem,pprime,q_target,model=None,psi0=None,*,
                        options=None,parallel=None,inner_options=None):
    """Direct FFprime inverse solve; every trial solves the full GS problem."""
    outer=_options(PQOptions,options); parallel=_options(ParallelOptions,parallel)
    inner=problem.picard if inner_options is None else _options(PicardOptions,inner_options)
    pp=_vector(pprime,"pprime",len(problem.psi_n))
    if model is None:
        model=initialize_ffprime_model(problem.psi_n,pp,q_target,r0=problem.r0,
                                      b0=problem.b0,a_minor=problem.a_minor)
    if not isinstance(model,FFprimeModel): raise TypeError("model must be FFprimeModel.")
    canonical=normalize_q_target(q_target,problem.psi_n)
    q_grid=np.r_[0.,problem.q_options.interior_psi_n,1.]
    target=evaluate_q_target(canonical,q_grid)
    guess=_vector(problem.initial_psi if psi0 is None else psi0,"psi0",problem.mesh.n_dof)
    scales=list(dict.fromkeys(outer.initial_scales)) if np.any(model.coefficients) else [1.]
    start=time.perf_counter()
    try:
        from threadpoolctl import threadpool_limits
    except ImportError as exc:
        raise ImportError("Install threadpoolctl to control BLAS threads.") from exc
    history=[]; last_update={}; evaluations=0; initial_trials=[]
    with threadpool_limits(limits=parallel.blas_threads):
        runtime=_Runtime(problem,pp,model,target,inner,outer)
        pool=_CandidatePool(runtime,parallel,max(model.n_basis,len(scales)))
        try:
            seeds=pool.map("candidate",[(scale*model.coefficients,guess,"",0) for scale in scales])
            evaluations+=len(seeds)
            initial_trials=[dict(scale=float(scale),valid=c.valid,merit=c.merit,failure=c.failure)
                            for scale,c in zip(scales,seeds)]
            valid=[c for c in seeds if c.valid]
            if not valid:
                details="; ".join(f"scale {s:g}: {c.failure}" for s,c in zip(scales,seeds))
                raise RuntimeError("No valid initial GS equilibrium. "+details)
            current=min(valid,key=lambda c:c.merit)
            current_scale=scales[next(i for i,c in enumerate(seeds) if c is current)]
            converged=False; reason="maximum_iterations"
            def record(candidate,step=0.,damping=np.nan,line_scale=np.nan):
                entry=dict(iteration=len(history)+1,merit=candidate.merit,
                           max_relative_q=candidate.metrics.max_abs_relative,
                           rms_log_q=candidate.metrics.rms_log,coefficients=candidate.coefficients.copy(),
                           coefficient_step=step,damping=damping,line_scale=line_scale,
                           inner_iterations=candidate.inner.iterations,evaluations=evaluations)
                history.append(entry)
                if outer.verbose:
                    print(f"q-solve {entry['iteration']:2d}: max={entry['max_relative_q']:.3e}, "
                          f"rms_log={entry['rms_log_q']:.3e}, candidate evaluations={evaluations}",flush=True)
            record(current,line_scale=float(current_scale))
            converged=_q_converged(current,outer) and len(history)>=outer.minimum_iterations
            if converged: reason="converged"
            while not converged and len(history)<outer.max_iterations:
                jobs=[(j,current.coefficients,current.state.psi,current.residual,
                       current.metrics.valid_mask,current.state.axis_mode,np.sign(current.state.dpsi))
                      for j in range(model.n_basis)]
                differences=pool.map("difference",jobs)
                evaluations+=sum(item[1]["evaluations"] for item in differences)
                if any(not item[1]["success"] for item in differences):
                    reason="FFprime_response_failure"
                    last_update=dict(finite_difference=tuple(item[1] for item in differences))
                    break
                response=np.column_stack([item[0] for item in differences])
                if not np.all(np.isfinite(response)):
                    reason="nonfinite_FFprime_response"; break
                damping=outer.initial_damping
                accepted=False; best=None; best_step=None; best_scale=np.nan; best_prediction=np.nan
                for _ in range(outer.maximum_damping_attempts):
                    try: step,prediction=_lm_step(current,response,damping,outer)
                    except np.linalg.LinAlgError:
                        damping*=outer.damping_growth; continue
                    size=np.linalg.norm(step,np.inf)
                    if size>outer.max_coefficient_step: step*=outer.max_coefficient_step/size
                    if np.linalg.norm(step,np.inf)<=outer.minimum_coefficient_step:
                        damping*=outer.damping_growth; continue
                    for line_index in range(outer.maximum_line_search_steps):
                        scale=.5**line_index
                        trial=runtime.candidate((current.coefficients+scale*step,current.state.psi,
                                                  current.state.axis_mode,np.sign(current.state.dpsi)))
                        evaluations+=1
                        if trial.valid and (best is None or trial.merit<best.merit):
                            best,best_step,best_scale,best_prediction=trial,scale*step,scale,scale*prediction
                        decrease=max(outer.minimum_absolute_merit_decrease,
                                     outer.armijo_fraction*scale*max(prediction,0.))
                        if trial.valid and (_q_converged(trial,outer) or trial.merit<=current.merit-decrease):
                            accepted=True
                            if _q_converged(trial,outer):
                                best,best_step,best_scale,best_prediction=trial,scale*step,scale,scale*prediction
                            break
                    if accepted: break
                    damping*=outer.damping_growth
                if not accepted and best is not None and best.merit<current.merit*(1-outer.minimum_relative_merit_decrease):
                    accepted=True
                if not accepted:
                    reason="no_acceptable_FFprime_step"; break
                previous=current.merit; current=best
                last_update=dict(coefficient_step=best_step,line_search_scale=best_scale,
                                 predicted_merit_reduction=best_prediction,
                                 actual_merit_reduction=previous-current.merit,response_matrix=response,
                                 finite_difference=tuple(item[1] for item in differences),damping=damping)
                record(current,float(np.linalg.norm(best_step,np.inf)),damping,best_scale)
                converged=_q_converged(current,outer) and len(history)>=outer.minimum_iterations
                if converged: reason="converged"
            parallel_info=dict(requested_workers=parallel.max_workers,workers=pool.workers,
                process_pool_started=pool.executor is not None,pids=tuple(sorted(pool.pids)),
                start_method=parallel.start_method,blas_threads=parallel.blas_threads,
                available_cpus=available_cpus(),rhs_cache_enabled=problem.rhs_cache is not None)
            def matrix_bytes(a):
                return a.data.nbytes+a.indices.nbytes+a.indptr.nbytes
            numeric_bytes=sum(matrix_bytes(a) for a in (problem.k,problem.kii,
                    runtime.cache.operator,problem.factor.L,problem.factor.U))
            if problem.rhs_cache is not None:
                numeric_bytes+=matrix_bytes(problem.rhs_cache.operator)
                numeric_bytes+=problem.rhs_cache.radius.nbytes+problem.rhs_cache.weight.nbytes
            numeric_bytes+=runtime.evaluator.inverse.nbytes+runtime.evaluator.origins.nbytes
            parallel_info["worker_numeric_array_mb_lower_bound"]=numeric_bytes/1024**2
        finally: pool.close()
    result=PQSolveResult(converged,reason,len(history),evaluations,time.perf_counter()-start,
        model.with_coefficients(current.coefficients),current.profile,current.inner,current.q_result,
        current.metrics,tuple(history),tuple(initial_trials),last_update,parallel_info)
    if not converged and outer.error_on_failure:
        raise RuntimeError(f"q-solve did not converge: {reason}, max relative q={current.metrics.max_abs_relative:.3e}.")
    return current.state,result


@dataclass(frozen=True)
class GSEquilibrium:
    mode: int
    problem: PreparedGSProblem
    state: GSState
    pprime: np.ndarray
    f_profile: object
    pressure: np.ndarray
    inner_result: PicardResult
    q_result: object = None
    pq_result: object = None

    def export_gfile(self, filename, **options):
        """Export a COCOS-1 gfile and provenance sidecar; see gs_fixed.geqdsk."""
        from .geqdsk import export_gfile
        return export_gfile(self, filename, **options)

    @property
    def mesh(self): return self.problem.mesh
    @property
    def psi(self): return self.state.psi
    @property
    def psi_n(self): return self.state.psi_n
    @property
    def axis_point(self):
        return self.q_result.axis.point if self.q_result is not None else self.state.axis_point
    @property
    def converged(self):
        return self.inner_result.converged and (self.pq_result is None or self.pq_result.converged)


@_single_blas_thread
def solve_gs(problem,pprime,*,ffprime=None,q_target=None,p_edge=0.,ffprime_model=None,
             initial=None,ffprime_options=None,outer_options=None,parallel=None,
             psi0=None,inner_options=None,compute_q=True):
    """Solve mode 1 (ffprime) or mode 3 (q_target) on a prepared problem."""
    if (ffprime is None)==(q_target is None):
        raise ValueError("Provide exactly one of ffprime or q_target.")
    pp=_vector(pprime,"pprime",len(problem.psi_n))
    if q_target is None:
        if ffprime_model is not None: raise ValueError("ffprime_model is only used in mode 3.")
        state,inner=solve_gs_inner(problem,pp,ffprime,psi0,inner_options)
        profile=reconstruct_g_from_ffprime(problem.psi_n,ffprime,problem.r0*problem.b0)
        if compute_q:
            evaluator=P3FieldEvaluator(problem.mesh)
            geometry=evaluate_q_geometry(problem.mesh,state,problem.q_options,evaluator=evaluator)
            q_result=evaluate_q_from_g(geometry,profile)
        else: q_result=None
        pq_result=None; mode=1
    else:
        model=ffprime_model
        if model is None:
            model=initialize_ffprime_model(problem.psi_n,pp,q_target,r0=problem.r0,b0=problem.b0,
                    a_minor=problem.a_minor,initial=initial,options=ffprime_options)
        state,pq_result=solve_pq_constraint(problem,pp,q_target,model,psi0,
                    options=outer_options,parallel=parallel,inner_options=inner_options)
        inner=pq_result.final_inner_result; profile=pq_result.final_f_profile
        q_result=pq_result.final_q_result; mode=3
    pressure=reconstruct_pressure(problem.psi_n,pp,p_edge)
    return GSEquilibrium(mode,problem,state,pp,profile,pressure,inner,q_result,pq_result)
