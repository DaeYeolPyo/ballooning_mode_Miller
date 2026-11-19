"""Solver and real process-pool tests: python -m test.test_solver."""
import os
import time
import numpy as np
from dataclasses import replace
from gs_fixed.solver import (
    PicardOptions,PQOptions,ParallelOptions,QOptions,prepare_gs_problem,
    solve_gs,solve_gs_inner,assemble_rhs,build_rhs_cache,available_cpus,
)
from gs_fixed.profiles import FFprimeModel,QTarget
from gs_fixed.post import P3FieldEvaluator,evaluate_q_geometry,build_q_contour_cache


def make_problem(nr=17,nz=21,contour_nr=96):
    theta=np.arange(80)*2*np.pi/80
    boundary=np.column_stack((3+.8*np.cos(theta+np.arcsin(.25)*np.sin(theta)),
                              1.2*np.sin(theta)))
    return prepare_gs_problem(boundary,r0=3,b0=3,a_minor=.8,nr=nr,nz=nz,n_profile=101,
        q_options=QOptions(contour_nr=contour_nr),
        picard=PicardOptions(max_iterations=150,update_tolerance=1e-9,residual_tolerance=1e-9))


def test_solver():
    problem=make_problem()
    s=problem.psi_n
    pp=-1e6*(1-s**1.5)**2
    c=np.array([2.,.2,-.1])
    from gs_fixed.profiles import evaluate_ffprime_model
    ff=evaluate_ffprime_model(FFprimeModel(s,c,9.)).ffprime
    direct=assemble_rhs(problem.mesh,problem.initial_psi,problem.quad,pp,ff,-.5)
    cached=problem.assemble_rhs(problem.initial_psi,pp,ff,-.5)
    np.testing.assert_allclose(cached,direct,rtol=3e-14,atol=1e-15)
    assert build_rhs_cache(problem.mesh,problem.quad,max_mb=.0001) is None

    t=time.perf_counter()
    reference=solve_gs(problem,pp,ffprime=ff)
    assert reference.converged
    assert reference.inner_result.update_history[-1]<1e-9
    assert reference.inner_result.residual_history[-1]<1e-9
    np.testing.assert_array_equal(reference.psi[problem.bc.boundary_nodes],0)
    np.testing.assert_array_equal(reference.psi_n[problem.bc.boundary_nodes],1)
    eq_seconds=time.perf_counter()-t

    # Cached and element-loop Picard paths solve the same nonlinear equation.
    loop=replace(problem,rhs_cache=None)
    uncached,_=solve_gs_inner(loop,pp,ff)
    np.testing.assert_allclose(uncached.psi,reference.psi,rtol=2e-12,atol=2e-13)
    # Sign-reversed flux/axis with the same normalized profiles.
    minimum=replace(problem,picard=replace(problem.picard,axis_mode="min"),
                     initial_psi=-problem.initial_psi)
    min_state,min_result=solve_gs_inner(minimum,pp,ff)
    assert min_result.converged and min_state.dpsi>0
    np.testing.assert_allclose(min_state.psi,-reference.psi,rtol=3e-12,atol=3e-13)

    # Cubic field values and gradients on the original constrained elements.
    evaluator=P3FieldEvaluator(problem.mesh)
    centers=problem.mesh.p1_points[problem.mesh.p1_elements].mean(axis=1)
    R,Z=problem.mesh.points.T
    for shift in (0.,3.):
        values=(R-shift)**3+(R-shift)*Z**2+2*Z
        val,grad,ids=evaluator.evaluate(values,centers)
        r,z=centers.T
        r=r-shift
        expected=np.column_stack((3*r*r+z*z,2*r*z+2))
        np.testing.assert_allclose(val,r**3+r*z**2+2*z,rtol=3e-14,atol=1e-13)
        # Derivatives amplify floating-point nodal errors on skinny triangles.
        bound=100*np.finfo(float).eps*np.linalg.norm(
            evaluator.inverse[ids],ord=np.inf,axis=(1,2))*np.max(
                abs(values[problem.mesh.elements[ids]]),axis=1)
        assert np.all(np.max(abs(grad-expected),axis=1)<=bound)
        if shift:
            np.testing.assert_allclose(grad,expected,rtol=3e-13,atol=3e-12)

    model=FFprimeModel(s,np.array([1.,0.,0.]),9.)
    target=QTarget(reference.q_result.psi_n,reference.q_result.q)
    options=PQOptions(max_iterations=12,q_relative_tolerance=1e-4,
        q_rms_log_tolerance=1e-4,regularization=0.,verbose=False)
    t=time.perf_counter()
    serial=solve_gs(problem,pp,q_target=target,ffprime_model=model,
                    outer_options=options,parallel=ParallelOptions(max_workers=1))
    serial_seconds=time.perf_counter()-t
    t=time.perf_counter()
    parallel=solve_gs(problem,pp,q_target=target,ffprime_model=model,
                    outer_options=options,parallel=ParallelOptions(max_workers=3))
    parallel_seconds=time.perf_counter()-t
    assert serial.converged and parallel.converged
    np.testing.assert_allclose(parallel.psi,serial.psi,rtol=2e-11,atol=2e-12)
    np.testing.assert_allclose(parallel.pq_result.final_model.coefficients,
                                serial.pq_result.final_model.coefficients,rtol=2e-10,atol=2e-10)
    assert serial.pq_result.function_evaluations==parallel.pq_result.function_evaluations
    assert serial.pq_result.outer_iterations==parallel.pq_result.outer_iterations
    if available_cpus()>=2:
        assert parallel.pq_result.parallel_info["process_pool_started"]
        assert any(pid!=os.getpid() for pid in parallel.pq_result.parallel_info["pids"])
    assert not serial.pq_result.parallel_info["process_pool_started"]

    # Target matches a zero-FFprime seed: lazy pool avoids spawning workers.
    zero=FFprimeModel(s,np.zeros(3),9.)
    exact=solve_gs(problem,pp,ffprime=np.zeros_like(s))
    matched=solve_gs(problem,pp,q_target=QTarget(exact.q_result.psi_n,exact.q_result.q),
        ffprime_model=zero,outer_options=options,parallel=ParallelOptions(max_workers=3))
    assert matched.converged and matched.pq_result.outer_iterations==1
    assert not matched.pq_result.parallel_info["process_pool_started"]

    # Explicit failure status and bad-input paths.
    from gs_fixed.solver import PicardOptions
    state,result=solve_gs_inner(problem,pp,ff,
        options=replace(problem.picard,max_iterations=1),emit_warning=False)
    assert not result.converged
    def fails(fn,*args,**kwargs):
        try: fn(*args,**kwargs)
        except (ValueError,TypeError,RuntimeError): return
        raise AssertionError("Invalid input accepted.")
    fails(solve_gs,problem,pp)
    fails(solve_gs,problem,pp,ffprime=ff,q_target=target)
    fails(solve_gs_inner,problem,pp,ff,options=replace(problem.picard,max_iterations=1,error_on_failure=True))
    fails(ParallelOptions,max_workers=-1)
    fails(PicardOptions,omega=1.1)
    fails(PQOptions,fit_range=(.9,.1))
    fails(solve_gs,problem,pp,q_target=-reference.q_result.q,ffprime_model=model,
          outer_options=options,parallel=ParallelOptions(max_workers=1))
    print("Solver, cache, field, branch, mode-1/mode-3 and process-equivalence checks passed.")
    print("DOFs",problem.mesh.n_dof,"mode1",eq_seconds,"serial q",serial_seconds,
           "parallel q",parallel_seconds,"seconds")
    print("q error",parallel.pq_result.final_metrics.max_abs_relative,
           "GS evaluations",parallel.pq_result.function_evaluations,
           "coefficients",parallel.pq_result.final_model.coefficients)
    print("Parallel diagnostics:",parallel.pq_result.parallel_info)


if __name__=="__main__":
    test_solver()
