"""P3 field evaluation and q geometry needed by the GS/q solvers.

Cartesian grids discover contours only. All gradients are evaluated from
the P3 field on its original affine triangles, never from grid differences.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from matplotlib.tri import Triangulation
from matplotlib.path import Path
import contourpy
from .fem import construct_p3_shape_functions
from .mesh import P3Mesh
from .profiles import validate_profile_grid, FProfile


@dataclass(frozen=True)
class QOptions:
    interior_psi_n: tuple = tuple(np.unique(np.r_[np.arange(1,17)*.05,
                                                   .825,.85,.875,.9,.925,.95]))
    contour_nr: int = 500
    contour_nz: int = 0

    def __post_init__(self):
        s=np.asarray(self.interior_psi_n,dtype=float)
        if s.ndim!=1 or len(s)<3 or not np.all(np.isfinite(s)) or np.any(np.diff(s)<=0) or np.any(s<=0) or np.any(s>=1):
            raise ValueError("At least three ordered q levels in (0,1) are required.")
        for name,minimum in (("contour_nr",16),("contour_nz",0)):
            value=getattr(self,name)
            if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or value<minimum:
                raise ValueError(f"Invalid {name}.")
        if 0<self.contour_nz<16:
            raise ValueError("contour_nz must be zero (automatic) or >=16.")
        object.__setattr__(self,"interior_psi_n",tuple(s))


@dataclass(frozen=True)
class AxisInfo:
    point: np.ndarray
    nodal_point: np.ndarray
    psi: float
    hessian: np.ndarray
    q_per_f: float
    axis_mode: str
    fit_residual: float
    patch_scale: float
    patch_node_count: int


@dataclass(frozen=True)
class QGeometry:
    psi_n: np.ndarray
    q_per_f: np.ndarray
    axis: AxisInfo
    minimum_grad_psi: np.ndarray
    contour_area: np.ndarray
    contour_point_count: np.ndarray
    grid_shape: tuple


@dataclass(frozen=True)
class QResult:
    psi_n: np.ndarray
    q: np.ndarray
    f: np.ndarray
    g: np.ndarray
    q_per_f: np.ndarray
    axis: AxisInfo
    geometry: QGeometry
    boundary_extrapolation_spread: float


@dataclass(frozen=True)
class QContourCache:
    operator: object
    inside: np.ndarray
    r: np.ndarray
    z: np.ndarray
    mesh_points: np.ndarray
    mesh_elements: np.ndarray
    solution_elements: np.ndarray


class P3FieldEvaluator:
    """A locator for the supplied constrained triangulation."""
    def __init__(self,mesh):
        if not isinstance(mesh,P3Mesh):
            raise TypeError("mesh must be P3Mesh.")
        self.mesh=mesh
        p,t=mesh.p1_points,mesh.p1_elements
        self.triangulation=Triangulation(p[:,0],p[:,1],t)
        self.finder=self.triangulation.get_trifinder()
        self.edge_max=np.zeros(len(mesh.edges),dtype=np.int64)
        self.vertex_max=np.zeros(len(p),dtype=np.int64)
        np.maximum.at(self.edge_max,mesh.element_edges.ravel(),np.repeat(np.arange(len(t)),3))
        np.maximum.at(self.vertex_max,t.ravel(),np.repeat(np.arange(len(t)),3))
        self.coordinate_scale=max(1.,float(np.max(np.abs(p))))
        vertices=p[t]
        self.origins=vertices[:,0]
        e1=vertices[:,1]-vertices[:,0]
        e2=vertices[:,2]-vertices[:,0]
        det=e1[:,0]*e2[:,1]-e2[:,0]*e1[:,1]
        if np.any(det<=0):
            raise ValueError("Affine elements must be counter-clockwise.")
        self.inverse=np.empty((len(t),2,2))
        self.inverse[:,0,0]=e2[:,1]/det
        self.inverse[:,0,1]=-e2[:,0]/det
        self.inverse[:,1,0]=-e1[:,1]/det
        self.inverse[:,1,1]=e1[:,0]/det

    def locate(self,points):
        x=np.asarray(points,dtype=float)
        if x.ndim!=2 or x.shape[1]!=2 or not np.all(np.isfinite(x)):
            raise ValueError("Query points must be finite (n,2).")
        ids=np.asarray(self.finder(x[:,0],x[:,1]),dtype=np.int64)
        inside=ids>=0
        local=np.empty((np.count_nonzero(inside),2))
        current=ids[inside]
        local[:]=np.einsum("nij,nj->ni",self.inverse[current],
                           x[inside]-self.origins[current])
        # MATLAB pointLocation selects the largest ID on a shared boundary.
        # P3 gradients are discontinuous there, so this tie rule matters.
        bary=np.column_stack((1-local.sum(axis=1),local))
        tolerance=64*np.finfo(float).eps*np.maximum(
            1.,np.linalg.norm(self.inverse[current],ord=np.inf,axis=(1,2))*self.coordinate_scale)
        on_edge=np.abs(bary)<=tolerance[:,None]
        chosen=current.copy()
        for opposite,local_edge in ((0,1),(1,2),(2,0)):
            maximum=self.edge_max[self.mesh.element_edges[current,local_edge]]
            chosen=np.where(on_edge[:,opposite],np.maximum(chosen,maximum),chosen)
        at_vertex=on_edge.sum(axis=1)>=2
        if np.any(at_vertex):
            local_vertex=np.argmax(bary[at_vertex],axis=1)
            vertex=self.mesh.p1_elements[current[at_vertex],local_vertex]
            chosen[at_vertex]=self.vertex_max[vertex]
        if np.any(chosen!=current):
            ids[inside]=chosen
            local[:]=np.einsum("nij,nj->ni",self.inverse[chosen],
                               x[inside]-self.origins[chosen])
        return ids,inside,local

    def evaluate(self,psi,points):
        psi=np.asarray(psi,dtype=float)
        if psi.shape!=(self.mesh.n_dof,) or not np.all(np.isfinite(psi)):
            raise ValueError("psi must contain one finite value per P3 node.")
        ids,inside,local=self.locate(points)
        value=np.full(len(ids),np.nan)
        gradient=np.full((len(ids),2),np.nan)
        if np.any(inside):
            n,dn=construct_p3_shape_functions(local[:,0],local[:,1])
            samples=psi[self.mesh.elements[ids[inside]]]
            value[inside]=np.einsum("ni,ni->n",n,samples)
            # A constant nodal offset has exactly zero gradient.
            reference_gradient=np.einsum("nki,ni->nk",dn,samples-samples[:,:1])
            gradient[inside]=np.einsum("nji,nj->ni",self.inverse[ids[inside]],
                                        reference_gradient)
        return value,gradient,ids

    def value_operator(self,points):
        ids,inside,local=self.locate(points)
        n,_=construct_p3_shape_functions(local[:,0],local[:,1])
        rows=np.repeat(np.flatnonzero(inside),10)
        columns=self.mesh.elements[ids[inside]].ravel()
        matrix=sparse.coo_matrix((n.ravel(),(rows,columns)),
                                 shape=(len(ids),self.mesh.n_dof)).tocsr()
        return matrix,inside


def eval_p3_field(mesh,psi,points,evaluator=None):
    evaluator=P3FieldEvaluator(mesh) if evaluator is None else evaluator
    if evaluator.mesh is not mesh:
        raise ValueError("Evaluator belongs to a different mesh.")
    return evaluator.evaluate(psi,points)


def build_q_contour_cache(mesh,options=None,evaluator=None):
    options=QOptions() if options is None else options
    evaluator=P3FieldEvaluator(mesh) if evaluator is None else evaluator
    if evaluator.mesh is not mesh:
        raise ValueError("Evaluator belongs to a different mesh.")
    p=mesh.p1_points
    r=np.linspace(p[:,0].min(),p[:,0].max(),options.contour_nr)
    nz=options.contour_nz or max(200,int(np.floor(options.contour_nr*np.ptp(p[:,1])/np.ptp(p[:,0])+.5)))
    z=np.linspace(p[:,1].min(),p[:,1].max(),nz)
    rr,zz=np.meshgrid(r,z)
    operator,inside=evaluator.value_operator(np.column_stack((rr.ravel(),zz.ravel())))
    return QContourCache(operator,inside,r,z,mesh.p1_points,mesh.p1_elements,mesh.elements)


def recover_axis(mesh,psi,axis_node,axis_mode,evaluator):
    center=mesh.points[axis_node]
    central=np.any(mesh.elements==axis_node,axis=1)
    nodes=np.unique(mesh.elements[central])
    expanded=np.any(np.isin(mesh.elements,nodes),axis=1)
    patch=np.unique(mesh.elements[expanded])
    displacement=mesh.points[patch]-center
    scale=np.max(np.linalg.norm(displacement,axis=1))
    if not np.isfinite(scale) or scale<=0:
        raise ValueError("Invalid axis recovery patch.")
    x,y=(displacement/scale).T
    a=np.column_stack((np.ones(len(x)),x,y,.5*x*x,x*y,.5*y*y))
    weights=1/(1+(np.hypot(x,y)/.5)**4)
    aw=a*np.sqrt(weights)[:,None]
    c,_,rank,_=np.linalg.lstsq(aw,psi[patch]*np.sqrt(weights),rcond=None)
    if rank<6:
        raise ValueError("Axis recovery patch has insufficient rank.")
    h=np.array([[c[3],c[4]],[c[4],c[5]]])/scale**2
    eigenvalues=np.linalg.eigvalsh(h)
    mode=axis_mode
    if mode=="auto":
        mode="max" if np.all(eigenvalues<0) else "min"
    if mode not in ("max","min") or (mode=="max" and np.any(eigenvalues>=0)) or (mode=="min" and np.any(eigenvalues<=0)):
        raise ValueError("Recovered axis Hessian has the wrong definiteness.")
    shift=-np.linalg.solve(h,c[1:3]/scale)
    point=center+shift
    if np.linalg.norm(shift)>.5*scale or point[0]<=0 or evaluator.finder(*point)<0:
        raise ValueError("Recovered magnetic axis lies outside its valid patch.")
    u,v=shift/scale
    axis_psi=float(np.array([1,u,v,.5*u*u,u*v,.5*v*v])@c)
    factor=1/(point[0]*np.sqrt(np.linalg.det(h)))
    denominator=max(np.linalg.norm(psi[patch]-np.mean(psi[patch])),
                    np.spacing(max(1.,np.linalg.norm(psi[patch]))))
    return AxisInfo(point,center.copy(),axis_psi,h,float(factor),mode,
                    float(np.linalg.norm(a@c-psi[patch])/denominator),float(scale),len(patch))


def evaluate_q_geometry(mesh,state,options=None,*,evaluator=None,cache=None):
    options=QOptions() if options is None else options
    evaluator=P3FieldEvaluator(mesh) if evaluator is None else evaluator
    if evaluator.mesh is not mesh:
        raise ValueError("Evaluator belongs to a different mesh.")
    cache=build_q_contour_cache(mesh,options,evaluator) if cache is None else cache
    if cache.operator.shape[1]!=mesh.n_dof or len(cache.r)!=options.contour_nr:
        raise ValueError("Incompatible contour cache.")
    for cached,actual in ((cache.mesh_points,mesh.p1_points),
                          (cache.mesh_elements,mesh.p1_elements),
                          (cache.solution_elements,mesh.elements)):
        if cached is not actual and not np.array_equal(cached,actual):
            raise ValueError("Contour cache belongs to a different mesh.")
    if options.contour_nz and len(cache.z)!=options.contour_nz:
        raise ValueError("Contour cache uses a different Z grid.")
    axis=recover_axis(mesh,state.psi,state.axis_node,state.axis_mode,evaluator)
    grid=np.asarray(cache.operator@state.psi)
    grid[~cache.inside]=np.nan
    grid=((grid-state.psi_axis)/state.dpsi).reshape(len(cache.z),len(cache.r))
    contours=contourpy.contour_generator(cache.r,cache.z,np.ma.masked_invalid(grid),
                                         name="serial",line_type="Separate",corner_mask=False)
    factors=[]; minimum=[]; areas=[]; counts=[]
    domain_scale=max(np.ptp(cache.r),np.ptp(cache.z),1.)
    for level in options.interior_psi_n:
        candidates=[]
        for curve in contours.lines(level):
            if len(curve)<4 or np.linalg.norm(curve[-1]-curve[0])>1e-8*domain_scale:
                continue
            curve=curve.copy()
            curve[-1]=curve[0]
            if Path(curve).contains_point(axis.point):
                shifted=curve-curve[0]
                area=.5*abs(np.sum(shifted[:-1,0]*shifted[1:,1]-shifted[1:,0]*shifted[:-1,1]))
                candidates.append((area,curve))
        if not candidates:
            raise ValueError(f"No closed axis-enclosing contour at psi_n={level:.6g}.")
        area,curve=max(candidates,key=lambda pair:pair[0])
        mid=.5*(curve[1:]+curve[:-1])
        dl=np.linalg.norm(np.diff(curve,axis=0),axis=1)
        _,gradient,_=evaluator.evaluate(state.psi,mid)
        magnitude=np.linalg.norm(gradient,axis=1)
        if np.any(~np.isfinite(magnitude)) or np.any(magnitude<=0) or np.any(mid[:,0]<=0):
            raise ValueError(f"Invalid P3 gradient at psi_n={level:.6g}.")
        integral=np.sum(dl/(mid[:,0]*magnitude))
        if not np.isfinite(integral) or integral<=0:
            raise ValueError("Invalid flux-surface integral.")
        factors.append(integral/(2*np.pi)); minimum.append(magnitude.min())
        areas.append(area); counts.append(len(curve))
    return QGeometry(np.r_[0.,options.interior_psi_n,1.],
                      np.r_[axis.q_per_f,factors,np.nan],axis,
                      np.asarray(minimum),np.asarray(areas),np.asarray(counts),grid.shape)


def evaluate_q_from_g(geometry,profile):
    if not isinstance(profile,FProfile):
        raise TypeError("profile must be FProfile.")
    s=validate_profile_grid(profile.psi_n)
    f=np.interp(geometry.psi_n,s,profile.f)
    interior=geometry.psi_n[1:-1]
    q_inside=f[1:-1]*geometry.q_per_f[1:-1]
    nfit=min(5,len(interior))
    x=interior[-nfit:]-1
    linear=np.polyfit(x,q_inside[-nfit:],1)[-1]
    edge=np.polyfit(x,q_inside[-nfit:],2)[-1] if nfit>=4 else linear
    spread=float(abs(edge-linear)) if nfit>=4 else np.nan
    q=np.r_[f[0]*geometry.q_per_f[0],q_inside,edge]
    if not np.all(np.isfinite(q)) or np.any(q==0):
        raise ValueError("Computed q must be finite and nonzero.")
    q_per_f=geometry.q_per_f.copy()
    q_per_f[-1]=edge/f[-1]
    return QResult(geometry.psi_n.copy(),q,f,f*f,q_per_f,geometry.axis,geometry,spread)

# Extensible global physical quantities, evaluated after the GS solve.

from scipy.interpolate import interp1d
from .fem import MU0, precompute_quadrature


@dataclass(frozen=True)
class PostData:
    """Shared integration samples for one equilibrium, in SI units.

    Arrays have shape (n_elements, n_quadrature_points). ``area_weights``
    represent dR*dZ; ``volume_weights`` represent 2*pi*R*dR*dZ. Weights may
    be negative for some valid Dunavant rules. Never replace them by abs().

    psi_n uses the SAME nodal axis and signed span as the GS source.
    Pressure/F use psi_n clipped to [0,1] within the requested tolerance;
    source derivatives use linear extrapolation, exactly as FEM assembly.
    Native field convention: BR=-psi_Z/R, BZ=psi_R/R, Bphi=F/R.
    """

    r: np.ndarray
    z: np.ndarray
    psi: np.ndarray
    psi_n: np.ndarray
    psi_n_profile: np.ndarray
    psi_r: np.ndarray
    psi_z: np.ndarray
    pressure: np.ndarray
    f: np.ndarray
    dp_dpsi: np.ndarray
    ffprime: np.ndarray
    j_phi: np.ndarray
    area_weights: np.ndarray
    volume_weights: np.ndarray
    r0: float
    b0: float
    a_minor: float
    r_axis: float
    z_axis: float
    b_axis: float
    psi_axis: float
    dpsi: float
    ip_from_rhs: float
    quadrature_degree: int
    converged: bool

    @property
    def b_r(self):
        return -self.psi_z / self.r

    @property
    def b_z(self):
        return self.psi_r / self.r

    @property
    def b_phi(self):
        return self.f / self.r

    @property
    def bp_squared(self):
        return (self.psi_r**2 + self.psi_z**2) / self.r**2

    @property
    def b_squared(self):
        return self.bp_squared + self.b_phi**2

    def _integrand(self, values):
        values = np.asarray(values)
        if not np.isrealobj(values):
            raise ValueError("Integrands must be real.")
        values = np.asarray(values, dtype=float)
        if values.ndim != 0 and values.shape != self.r.shape:
            raise ValueError("Use a scalar or an array with shape data.r.shape.")
        if not np.all(np.isfinite(values)):
            raise ValueError("Integrands must be finite at every quadrature point.")
        return values

    def integrate_area(self, values=1.0):
        """Integral over the poloidal cross-section: integral values dR dZ."""
        return float(np.sum(self._integrand(values) * self.area_weights))

    def integrate_volume(self, values=1.0):
        """Axisymmetric volume integral: integral values * 2*pi*R dR dZ."""
        return float(np.sum(self._integrand(values) * self.volume_weights))

    def volume_average(self, values):
        """Return integral values dV / integral dV."""
        volume = self.integrate_volume()
        if not np.isfinite(volume) or volume <= 0:
            raise ValueError("Plasma volume must be finite and positive.")
        return self.integrate_volume(values) / volume


def _post_scalar(value, name, *, positive=False, nonzero=False):
    raw = np.asarray(value)
    if raw.ndim != 0 or not np.isrealobj(raw) or raw.dtype.kind not in "iuf":
        raise ValueError(name + " must be a real scalar.")
    value = float(raw)
    if not np.isfinite(value) or (positive and value <= 0) or (nonzero and value == 0):
        raise ValueError(name + " has an invalid value.")
    return value


def _post_vector(value, name, size):
    raw = np.asarray(value)
    if not np.isrealobj(raw):
        raise ValueError(name + " must be real.")
    value = np.asarray(raw, dtype=float)
    if value.shape != (size,) or not np.all(np.isfinite(value)):
        raise ValueError(name + " must be a finite vector of the required length.")
    return value


def _post_inputs(equilibrium, allow_unconverged):
    # Local import avoids a post -> solver -> post import cycle.
    from .solver import GSEquilibrium
    if not isinstance(equilibrium, GSEquilibrium):
        raise TypeError("equilibrium must be a GSEquilibrium returned by solve_gs().")
    if not isinstance(allow_unconverged, (bool, np.bool_)):
        raise ValueError("allow_unconverged must be boolean.")
    if not equilibrium.converged and not allow_unconverged:
        raise ValueError("Equilibrium did not converge; set allow_unconverged=True for diagnostics.")
    problem, state = equilibrium.problem, equilibrium.state
    grid = validate_profile_grid(problem.psi_n, require_uniform=True)
    n = len(grid)
    psi = _post_vector(state.psi, "psi", problem.mesh.n_dof)
    axis = _post_scalar(state.psi_axis, "psi_axis")
    span = _post_scalar(state.dpsi, "dpsi", nonzero=True)
    pressure = _post_vector(equilibrium.pressure, "pressure", n)
    if np.min(pressure) < -1.e-10 * max(1., float(np.max(abs(pressure)))):
        raise ValueError("Pressure is negative; check p_edge and dp/dpsi_N.")
    pp = _post_vector(equilibrium.pprime, "pprime", n)
    profile = equilibrium.f_profile
    fgrid = validate_profile_grid(profile.psi_n)
    if fgrid.shape != grid.shape or not np.allclose(fgrid, grid, rtol=0., atol=1.e-13):
        raise ValueError("F profile must use the common problem.psi_n grid.")
    f = _post_vector(profile.f, "F", n)
    ff = _post_vector(profile.ffprime, "FFprime", n)
    if np.any(f == 0) or np.any(np.sign(f) != np.sign(f[0])):
        raise ValueError("F must remain nonzero and keep its sign.")
    # The derivatives are dp/dpsi_N and F*dF/dpsi_N at input.
    derivatives = interp1d(grid, np.column_stack((pp, ff)) / span, axis=0,
                          bounds_error=False, fill_value="extrapolate", assume_sorted=True)
    return grid, psi, axis, span, pressure, f, derivatives


def _post_profile_samples(psi_n, grid, pressure, f, derivatives, tolerance):
    finite = np.isfinite(psi_n)
    if np.any(psi_n[finite] < -tolerance) or np.any(psi_n[finite] > 1 + tolerance):
        raise ValueError("psi_N is outside [0,1] beyond flux_range_tolerance.")
    clipped = np.clip(psi_n, 0., 1.)
    return clipped, np.interp(clipped, grid, pressure), np.interp(clipped, grid, f), derivatives(psi_n)


def build_post_data(equilibrium, *, quadrature_degree=None,
                    flux_range_tolerance=1.e-3, allow_unconverged=False,
                    element_chunk_size=256):
    """Sample a solved equilibrium ONCE for any number of physical quantities.

    The default rule is the solver's rule. A different degree is useful for
    integration convergence checks; then current from these samples need
    not equal current from the original solver RHS exactly.
    """
    grid, psi, axis, span, pressure, f, derivatives = _post_inputs(equilibrium, allow_unconverged)
    tolerance = _post_scalar(flux_range_tolerance, "flux_range_tolerance")
    if tolerance < 0:
        raise ValueError("flux_range_tolerance must be nonnegative.")
    if (isinstance(element_chunk_size, (bool, np.bool_)) or
            not isinstance(element_chunk_size, (int, np.integer)) or element_chunk_size < 1):
        raise ValueError("element_chunk_size must be a positive integer.")
    problem = equilibrium.problem
    mesh = problem.mesh
    quad = problem.quad if quadrature_degree is None else precompute_quadrature(quadrature_degree)
    shape = (len(mesh.elements), quad.nq)
    # Chunk temporary basis-gradient arrays; retained fields have shape (nt,nq).
    r, z, psi_q, psi_r, psi_z, darea = [np.empty(shape) for _ in range(6)]
    for first in range(0, shape[0], int(element_chunk_size)):
        section = slice(first, first + int(element_chunk_size))
        vertices = mesh.p1_points[mesh.p1_elements[section]]
        origin = vertices[:, 0]
        e1, e2 = vertices[:, 1] - origin, vertices[:, 2] - origin
        det = e1[:, 0]*e2[:, 1] - e2[:, 0]*e1[:, 1]
        if np.any(det <= 0):
            raise ValueError("Elements must have positive affine Jacobians.")
        r[section] = origin[:, 0, None] + e1[:, 0, None]*quad.xi + e2[:, 0, None]*quad.eta
        z[section] = origin[:, 1, None] + e1[:, 1, None]*quad.xi + e2[:, 1, None]*quad.eta
        local = psi[mesh.elements[section]]
        psi_q[section] = local @ quad.N.T
        # Subtracting a constant preserves gradients while reducing cancellation.
        shifted = local - local[:, :1]
        gxi, geta = shifted @ quad.dNdxi.T, shifted @ quad.dNdeta.T
        psi_r[section] = (e2[:, 1, None]*gxi - e1[:, 1, None]*geta) / det[:, None]
        psi_z[section] = (-e2[:, 0, None]*gxi + e1[:, 0, None]*geta) / det[:, None]
        darea[section] = det[:, None]*quad.w
    if np.any(r <= 0):
        raise ValueError("Quadrature radii must be positive.")
    psi_n = (psi_q - axis) / span
    clipped, p_q, f_q, source = _post_profile_samples(psi_n, grid, pressure, f, derivatives, tolerance)
    j_phi = r*source[..., 0] + source[..., 1]/(MU0*r)
    dvolume = 2*np.pi*r*darea
    if not np.isfinite(np.sum(dvolume)) or np.sum(dvolume) <= 0:
        raise ValueError("Plasma volume must be finite and positive.")
    axis_point = np.asarray(equilibrium.axis_point, dtype=float)
    if axis_point.shape != (2,) or not np.all(np.isfinite(axis_point)) or axis_point[0] <= 0:
        raise ValueError("Invalid magnetic-axis coordinates.")
    r0 = _post_scalar(problem.r0, "r0", positive=True)
    b0 = _post_scalar(problem.b0, "b0", nonzero=True)
    minor = _post_scalar(problem.a_minor, "a_minor", positive=True)
    rhs = _post_vector(equilibrium.state.final_rhs, "final_rhs", mesh.n_dof)
    arrays = (r, z, psi_q, psi_n, clipped, psi_r, psi_z, p_q, f_q,
              source[..., 0], source[..., 1], j_phi, darea, dvolume)
    for array in arrays:
        if not np.all(np.isfinite(array)):
            raise ValueError("Nonfinite postprocessing samples.")
        array.setflags(write=False)
    return PostData(*arrays, r0, b0, minor, float(axis_point[0]), float(axis_point[1]),
                    float(f[0]/axis_point[0]), axis, span, float(np.sum(rhs)/MU0),
                    int(quad.degree), bool(equilibrium.converged))


def evaluate_magnetic_field(equilibrium, r, z, *, evaluator=None,
                            flux_range_tolerance=1.e-3, allow_unconverged=False):
    """Return a dict of native P3 fields/profiles at broadcast-compatible R,Z.

    Entries retain the broadcast query shape. Outside the original FEM
    domain ``inside`` is False and numeric field entries are NaN.
    Reuse a P3FieldEvaluator for repeated spatial queries.
    """
    grid, psi, axis, span, pressure, f, derivatives = _post_inputs(equilibrium, allow_unconverged)
    tolerance = _post_scalar(flux_range_tolerance, "flux_range_tolerance")
    if tolerance < 0:
        raise ValueError("flux_range_tolerance must be nonnegative.")
    if not np.isrealobj(r) or not np.isrealobj(z):
        raise ValueError("Query coordinates must be real.")
    r, z = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(z, dtype=float))
    shape = r.shape
    value, gradient, ids = eval_p3_field(equilibrium.mesh, psi,
                                        np.column_stack((r.ravel(), z.ravel())), evaluator)
    psi_n = (value - axis) / span
    _, p, f_q, source = _post_profile_samples(psi_n, grid, pressure, f, derivatives, tolerance)
    inside = ids >= 0
    safe_r = np.where(inside, r.ravel(), np.nan)
    if np.any(safe_r[inside] <= 0):
        raise ValueError("Inside query radii must be positive.")
    br, bz, bphi = -gradient[:, 1]/safe_r, gradient[:, 0]/safe_r, f_q/safe_r
    bp = np.hypot(br, bz)
    result = dict(psi=value, psi_n=psi_n, psi_r=gradient[:, 0], psi_z=gradient[:, 1],
                  b_r=br, b_z=bz, b_phi=bphi, bp=bp, b=np.hypot(bp, bphi),
                  pressure=p, f=f_q, dp_dpsi=source[:, 0], ffprime=source[:, 1],
                  j_phi=safe_r*source[:, 0] + source[:, 1]/(MU0*safe_r), inside=inside)
    return {name: values.reshape(shape) for name, values in result.items()}


def compute_plasma_current(data):
    """Signed Ip [A] = integral_Omega (R*dp/dpsi + FFprime/(mu0*R)) dA."""
    return data.integrate_area(data.j_phi)


def compute_plasma_volume(data):
    """Plasma volume [m^3] = integral_Omega 2*pi*R dA."""
    return data.integrate_volume()


def compute_pressure_integral(data):
    """Return integral p dV [J]; this is NOT yet a thermal-energy definition."""
    return data.integrate_volume(data.pressure)


def compute_poloidal_magnetic_energy(data):
    """Wp [J] = integral Bp^2 dV / (2*mu0)."""
    return data.integrate_volume(data.bp_squared) / (2*MU0)


def _post_current(data, ip):
    value = compute_plasma_current(data) if ip is None else _post_scalar(ip, "ip")
    if not np.isfinite(value) or abs(value) <= np.finfo(float).eps*max(1., abs(value)):
        raise ValueError("A nonzero plasma current is required for beta_pol, beta_N and li3.")
    return value


def compute_beta_tor(data):
    """Dimensionless beta_tor = 2*mu0*<p>_V / B0^2 (reference B0)."""
    return 2*MU0*data.volume_average(data.pressure) / data.b0**2


def compute_beta_pol(data, *, ip=None):
    """MATLAB convention: beta_pol = 4*integral p dV / (R0*mu0*Ip^2)."""
    current = _post_current(data, ip)
    return 4*compute_pressure_integral(data) / (data.r0*MU0*current**2)


def compute_li3(data, *, ip=None):
    """li(3) = 2*integral Bp^2 dV / (R0*mu0^2*Ip^2)."""
    current = _post_current(data, ip)
    return 2*data.integrate_volume(data.bp_squared) / (data.r0*MU0**2*current**2)


def compute_beta_n(data, *, ip=None, beta_tor=None):
    """Conventional numeric beta_N = beta_tor[%]*a[m]*abs(B0[T])/abs(Ip[MA])."""
    current = _post_current(data, ip)
    beta = compute_beta_tor(data) if beta_tor is None else _post_scalar(beta_tor, "beta_tor")
    return (100*beta)*data.a_minor*abs(data.b0) / (abs(current)/1.e6)


def compute_global_parameters(data):
    """Collect current quantities into an extensible plain dict.

    Add new compute_xxx(data) functions without editing solver.py or
    PostData. Call them directly or add their returned values to this dict.
    beta_tor/beta_pol are fractions; beta_N uses the percent/MA convention.
    """
    if not isinstance(data, PostData):
        raise TypeError("Build PostData with build_post_data(equilibrium) first.")
    ip = compute_plasma_current(data)
    beta = compute_beta_tor(data)
    scale = max(abs(ip), abs(data.ip_from_rhs), np.finfo(float).tiny)
    return dict(ip=ip, volume=compute_plasma_volume(data),
                pressure_integral=compute_pressure_integral(data),
                pressure_volume_average=data.volume_average(data.pressure),
                bp2_volume_integral=data.integrate_volume(data.bp_squared),
                poloidal_magnetic_energy=compute_poloidal_magnetic_energy(data),
                beta_tor=beta, beta_pol=compute_beta_pol(data, ip=ip),
                beta_n=compute_beta_n(data, ip=ip, beta_tor=beta), li3=compute_li3(data, ip=ip),
                r_axis=data.r_axis, z_axis=data.z_axis, b_axis=data.b_axis,
                ip_from_rhs=data.ip_from_rhs, ip_rhs_relative_difference=abs(ip-data.ip_from_rhs)/scale,
                quadrature_degree=data.quadrature_degree, converged=data.converged,
                psi_n_min=float(np.min(data.psi_n)), psi_n_max=float(np.max(data.psi_n)))

def plot_gs_equilibrium(eq, *, nr=240, q_target=None, quantities=None):
    """Create a figure; the caller chooses show/save/close."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from .profiles import normalize_q_target, evaluate_q_target

    cache = build_q_contour_cache(eq.mesh, QOptions(contour_nr=nr))
    values = (cache.operator @ eq.psi - eq.state.psi_axis) / eq.state.dpsi
    psi_n = np.ma.array(
        values.reshape(len(cache.z), len(cache.r)),
        mask=(~cache.inside | ~np.isfinite(values)).reshape(
            len(cache.z), len(cache.r)),
    )

    fig = plt.figure(figsize=(13, 7.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, width_ratios=(1.2, 1, 1))
    axes = {
        "flux": fig.add_subplot(grid[:, 0]),
        "pressure": fig.add_subplot(grid[0, 1]),
        "f": fig.add_subplot(grid[0, 2]),
        "q": fig.add_subplot(grid[1, 1]),
        "convergence": fig.add_subplot(grid[1, 2]),
    }

    ax = axes["flux"]
    filled = ax.contourf(
        cache.r, cache.z, psi_n, levels=np.linspace(0, 1, 21),
        cmap="viridis", extend="both", corner_mask=False,
    )
    ax.contour(
        cache.r, cache.z, psi_n, levels=np.arange(.1, 1, .1),
        colors="k", linewidths=.5, alpha=.5, corner_mask=False,
    )
    segments = eq.mesh.p1_points[eq.mesh.boundary_edges]
    ax.add_collection(LineCollection(segments, colors="k", linewidths=1.3))
    ax.plot(*eq.axis_point, marker="*", color="red", markersize=11,
            markeredgecolor="white", label="Magnetic axis")
    ax.set(xlabel="R [m]", ylabel="Z [m]", title="Flux surfaces")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(loc="best", fontsize=8)
    fig.colorbar(filled, ax=ax, label=r"$\psi_N$", shrink=.8)

    s = eq.problem.psi_n
    axes["pressure"].plot(s, eq.pressure / 1e3)
    axes["pressure"].set(title="Pressure", ylabel="p [kPa]")
    axes["f"].plot(eq.f_profile.psi_n, eq.f_profile.f)
    axes["f"].set(title="Toroidal-field function", ylabel="F [T m]")

    ax = axes["q"]
    if eq.q_result is not None:
        sq = eq.q_result.psi_n
        ax.plot(sq, eq.q_result.q, label="Computed")
        if q_target is not None:
            target = normalize_q_target(q_target, s)
            ax.plot(sq, evaluate_q_target(target, sq), "--", label="Target")
        ax.legend(fontsize=8)
    else:
        ax.text(.5, .5, "q was not computed", ha="center",
                va="center", transform=ax.transAxes)
    ax.set(title="Safety factor", ylabel="q")

    for key in ("pressure", "f", "q"):
        axes[key].set(xlabel=r"$\psi_N$", xlim=(0, 1))
        axes[key].grid(True, alpha=.25)

    ax = axes["convergence"]

    def log_curve(values, label):
        values = np.asarray(values, dtype=float)
        ax.semilogy(
            np.arange(1, values.size + 1),
            np.maximum(values, np.finfo(float).tiny), "o-",
            markersize=3, label=label,
        )

    if eq.pq_result is not None:
        history = eq.pq_result.history
        log_curve([h["max_relative_q"] for h in history], "Max relative q error")
        log_curve([h["rms_log_q"] for h in history], "RMS log q error")
        ax.set(title="Outer q-solve", xlabel="Outer iteration")
    else:
        log_curve(eq.inner_result.update_history, "Relative update")
        log_curve(eq.inner_result.residual_history, "Relative residual")
        ax.set(title="Inner GS solve", xlabel="Picard iteration")
    ax.set_ylabel("Error")
    ax.grid(True, which="both", alpha=.25)
    ax.legend(fontsize=8)

    status = "converged" if eq.converged else "NOT converged"
    title = f"GS equilibrium | mode {eq.mode} | {status}"
    if quantities is not None:
        title += (
            f"\nIp = {quantities['ip']/1e6:.3f} MA"
            f"   beta_tor = {100*quantities['beta_tor']:.2f}%"
            f"   beta_pol = {quantities['beta_pol']:.3f}"
            f"   beta_N = {quantities['beta_n']:.3f}"
            f"   li(3) = {quantities['li3']:.3f}"
        )
    fig.suptitle(title, fontsize=12)
    return fig, axes

# Append this block ONCE to gs_fixed/post.py. Existing post imports are used.
# Definition: S = -h dot curl(h), h = (B cross grad(psi))/|grad(psi)|^2.


class P3HessianEvaluator:
    """Native P3 values/gradients and analytical element Hessians.

    hessian_method='p3' returns the piecewise Hessian of the original field.
    'l2' projects ONLY the Hessian onto continuous P1 nodes for a recovered
    diagnostic. The latter is an approximation, not an exact element curl.
    Both methods retain the original constrained mesh and P3 gradient.
    """
    def __init__(self, mesh, psi, *, hessian_method="p3"):
        from .fem import P3_REFERENCE_NODES
        if hessian_method not in ("p3", "l2"):
            raise ValueError("hessian_method must be 'p3' or 'l2'.")
        self.mesh = mesh
        if not np.isrealobj(psi):
            raise ValueError("psi must be real.")
        self.psi = np.array(psi, dtype=float, copy=True)
        if self.psi.shape != (mesh.n_dof,) or not np.all(np.isfinite(self.psi)):
            raise ValueError("psi must be a finite vector with one value per P3 node.")
        self.psi.setflags(write=False)
        self.hessian_method = hessian_method
        self.locator = P3FieldEvaluator(mesh)
        x, y = P3_REFERENCE_NODES.T
        vandermonde = np.column_stack((np.ones(10), x, y, x*x, x*y, y*y,
                                      x**3, x*x*y, x*y*y, y**3))
        local = self.psi[mesh.elements]
        self.coefficients = np.linalg.solve(vandermonde, (local-local[:, :1]).T).T
        self.recovered_hessian = self._recover_hessian() if hessian_method == "l2" else None

    def _element_hessian(self, ids, local):
        c = self.coefficients[ids]
        x, y = local.T
        hxx = 2*c[:, 3] + 6*x*c[:, 6] + 2*y*c[:, 7]
        hxy = c[:, 4] + 2*x*c[:, 7] + 2*y*c[:, 8]
        hyy = 2*c[:, 5] + 2*x*c[:, 8] + 6*y*c[:, 9]
        inv = self.locator.inverse[ids]
        a, b, c, d = inv[:, 0, 0], inv[:, 0, 1], inv[:, 1, 0], inv[:, 1, 1]
        out = np.empty((len(ids), 2, 2))
        out[:, 0, 0] = a*a*hxx + 2*a*c*hxy + c*c*hyy
        out[:, 0, 1] = out[:, 1, 0] = a*b*hxx + (a*d+b*c)*hxy + c*d*hyy
        out[:, 1, 1] = b*b*hxx + 2*b*d*hxy + d*d*hyy
        return out

    def _recover_hessian(self):
        from scipy.sparse.linalg import spsolve
        t = self.mesh.p1_elements
        ids = np.repeat(np.arange(len(t)), 3)
        local = np.tile(np.array([[0., 0.], [1., 0.], [0., 1.]]), (len(t), 1))
        h = self._element_hessian(ids, local)
        values = np.column_stack((h[:, 0, 0], h[:, 0, 1], h[:, 1, 1])).reshape(-1, 3, 3)
        vertices = self.mesh.p1_points[t]
        e1, e2 = vertices[:, 1]-vertices[:, 0], vertices[:, 2]-vertices[:, 0]
        det = e1[:, 0]*e2[:, 1] - e1[:, 1]*e2[:, 0]
        mass = det[:, None, None]*(np.ones((3, 3))+np.eye(3))[None, :, :]/24
        rows, cols = np.repeat(t, 3, axis=1).ravel(), np.tile(t, (1, 3)).ravel()
        nv = len(self.mesh.p1_points)
        matrix = sparse.coo_matrix((mass.ravel(), (rows, cols)), shape=(nv, nv)).tocsc()
        rhs = np.zeros((nv, 3))
        local_rhs = np.einsum("eij,ejk->eik", mass, values)
        np.add.at(rhs, t.ravel(), local_rhs.reshape(-1, 3))
        recovered = spsolve(matrix, rhs)
        if not np.all(np.isfinite(recovered)):
            raise ValueError("Hessian recovery failed.")
        return recovered

    def evaluate(self, points):
        ids, inside, local = self.locator.locate(points)
        value = np.full(len(ids), np.nan)
        gradient = np.full((len(ids), 2), np.nan)
        hessian = np.full((len(ids), 2, 2), np.nan)
        if np.any(inside):
            element_ids = ids[inside]
            n, dn = construct_p3_shape_functions(local[:, 0], local[:, 1])
            samples = self.psi[self.mesh.elements[element_ids]]
            value[inside] = np.einsum("ni,ni->n", n, samples)
            reference_gradient = np.einsum("nki,ni->nk", dn, samples-samples[:, :1])
            gradient[inside] = np.einsum("nji,nj->ni", self.locator.inverse[element_ids], reference_gradient)
            if self.hessian_method == "p3":
                hessian[inside] = self._element_hessian(element_ids, local)
            else:
                bary = np.column_stack((1-local.sum(axis=1), local))
                nodal = self.recovered_hessian[self.mesh.p1_elements[element_ids]]
                h = np.einsum("ni,nij->nj", bary, nodal)
                hessian[inside, 0, 0] = h[:, 0]
                hessian[inside, 0, 1] = hessian[inside, 1, 0] = h[:, 1]
                hessian[inside, 1, 1] = h[:, 2]
        return value, gradient, hessian, ids


def _local_shear_formula(r, gradient, hessian, f, fprime, poloidal_sign):
    """Cylindrical expansion of the supplied NEGATIVE h.dot(curl(h))."""
    a, b = gradient[:, 0], gradient[:, 1]
    rr, rz, zz = hessian[:, 0, 0], hessian[:, 0, 1], hessian[:, 1, 1]
    g = a*a+b*b
    normal_second = (a*a*rr+2*a*b*rz+b*b*zz)/g
    return -poloidal_sign/r**2 * (fprime + f/g*(rr+zz-a/r-2*normal_second))


def _local_shear_f_profile(profile, psi_n, signed_span):
    """Integrate piecewise-linear F*dF/du exactly BETWEEN profile knots.

    This preserves the stored cumulative-trapezoid G values AND pairs F
    with its actual derivative. Differentiating an unrelated F interpolant
    while using FFprime from another interpolant would break this identity.
    """
    grid = validate_profile_grid(profile.psi_n)
    g = np.asarray(profile.g, dtype=float)
    ff = np.asarray(profile.ffprime, dtype=float)
    f_nodes = np.asarray(profile.f, dtype=float)
    if any(a.shape != grid.shape or not np.all(np.isfinite(a)) for a in (g, ff, f_nodes)):
        raise ValueError("Invalid F profile arrays.")
    if np.any(g <= 0) or np.any(f_nodes == 0) or np.any(np.sign(f_nodes) != np.sign(f_nodes[0])):
        raise ValueError("F must remain nonzero with a fixed sign.")
    step = np.diff(grid)
    error = np.diff(g)-(ff[1:]+ff[:-1])*step
    g_scale = max(float(np.max(abs(g))), np.finfo(float).tiny)
    if np.max(abs(g-f_nodes**2)) > 1e-10*g_scale:
        raise ValueError("Stored F and F^2 are inconsistent.")
    if np.max(abs(error)) > 1e-10*g_scale:
        raise ValueError("F^2 and FFprime are inconsistent with cumulative trapezoidal reconstruction.")
    i = np.clip(np.searchsorted(grid, psi_n, side="right")-1, 0, len(grid)-2)
    offset = psi_n-grid[i]
    slope = np.diff(ff)[i]/step[i]
    g_query = g[i]+2*ff[i]*offset+slope*offset**2
    if np.any(g_query <= 0) or not np.all(np.isfinite(g_query)):
        raise ValueError("F^2 becomes nonpositive between profile knots; refine/check the profile.")
    f = np.sign(f_nodes[0])*np.sqrt(g_query)
    fprime = (ff[i]+slope*offset)/(f*signed_span)
    return f, fprime


def compute_local_shear(equilibrium, r, z, *, hessian_method="p3",
                        poloidal_sign=1, psi_n_range=(.02, .98),
                        relative_gradient_floor=1e-8, evaluator=None):
    """Evaluate S [m^-3] from the user's first definition.

    B_R=-poloidal_sign*psi_Z/R, B_Z=poloidal_sign*psi_R/R, Bphi=F/R;
    psi is the NATIVE dimensional GS flux. Default +1 matches post.py's
    native magnetic-field convention. Reversing ONLY the poloidal B sign
    reverses S but preserves its zero set. This is not magnetic shear shat.

    p3: exact polynomial element Hessian (discontinuous across C0 elements).
    l2: recovered P1 Hessian for a separate convergence diagnostic; it is
        not an exact curl of the original piecewise field at finite mesh.
    """
    from .solver import GSEquilibrium
    if not isinstance(equilibrium, (GSEquilibrium, LocalShearReference)):
        raise TypeError("Use a solve_gs() result or load_local_shear_reference().")
    if not equilibrium.converged:
        raise ValueError("Refuse an unconverged GS equilibrium.")
    if np.ndim(poloidal_sign) != 0 or isinstance(poloidal_sign, (bool, np.bool_)) or poloidal_sign not in (-1, 1):
        raise ValueError("poloidal_sign must be +1 or -1.")
    limits = np.asarray(psi_n_range, dtype=float)
    if limits.shape != (2,) or not np.all(np.isfinite(limits)) or not 0 <= limits[0] < limits[1] <= 1:
        raise ValueError("psi_n_range must satisfy 0 <= lower < upper <= 1.")
    floor = float(relative_gradient_floor)
    span, axis = float(equilibrium.state.dpsi), float(equilibrium.state.psi_axis)
    minor = float(equilibrium.problem.a_minor)
    if not np.all(np.isfinite([floor, span, axis, minor])) or floor < 0 or span == 0 or minor <= 0:
        raise ValueError("Invalid gradient floor, flux span, axis or reference length.")
    if not np.isrealobj(r) or not np.isrealobj(z):
        raise ValueError("Query coordinates must be real.")
    r, z = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(z, dtype=float))
    shape = r.shape
    if evaluator is None:
        evaluator = P3HessianEvaluator(equilibrium.mesh, equilibrium.psi, hessian_method=hessian_method)
    elif (evaluator.mesh is not equilibrium.mesh or evaluator.hessian_method != hessian_method
          or not np.array_equal(evaluator.psi, equilibrium.psi)):
        raise ValueError("Evaluator must match the mesh, solution and Hessian method.")
    psi, gradient, hessian, ids = evaluator.evaluate(np.column_stack((r.ravel(), z.ravel())))
    psi_n = (psi-axis)/span
    magnitude = np.hypot(gradient[:, 0], gradient[:, 1])
    threshold = max(floor*abs(span)/minor, np.sqrt(np.finfo(float).tiny))
    valid = (ids >= 0) & (r.ravel() > 0) & np.isfinite(magnitude) & (magnitude > threshold)
    valid &= (psi_n >= limits[0]) & (psi_n <= limits[1])
    valid &= np.all(np.isfinite(hessian), axis=(1, 2))
    values = np.full(psi.shape, np.nan)
    f, fp = np.full(psi.shape, np.nan), np.full(psi.shape, np.nan)
    if np.any(valid):
        f[valid], fp[valid] = _local_shear_f_profile(equilibrium.f_profile, psi_n[valid], span)
        values[valid] = _local_shear_formula(r.ravel()[valid], gradient[valid],
                                           hessian[valid], f[valid], fp[valid], poloidal_sign)
        if not np.all(np.isfinite(values[valid])):
            raise ValueError("Nonfinite local shear; increase the axis/gradient exclusion.")
    out = dict(S=values, valid=valid, psi_n=psi_n, grad_psi=magnitude,
               f=f, dF_dpsi=fp, psi_r=gradient[:, 0], psi_z=gradient[:, 1],
               psi_rr=hessian[:, 0, 0], psi_rz=hessian[:, 0, 1], psi_zz=hessian[:, 1, 1])
    return {name:a.reshape(shape) for name, a in out.items()}


@dataclass(frozen=True)
class LocalShearMap:
    r: np.ndarray
    z: np.ndarray
    S: np.ndarray
    psi_n: np.ndarray
    valid: np.ndarray
    hessian_method: str
    poloidal_sign: int
    psi_n_range: tuple
    relative_gradient_floor: float
    units: str = "m^-3"


@dataclass(frozen=True)
class LocalShearReference:
    """Minimal read-only postprocessing view; not a restartable GS solve."""
    mesh: object
    psi: np.ndarray
    f_profile: object
    state: object
    problem: object
    axis_point: np.ndarray
    converged: bool


def load_local_shear_reference(reference_file, *, metadata_file=None):
    """Read reference_gs.npz + metadata.json saved by run_gs_miller_salpha.

    No global GS solve or Miller fit is rerun. Convergence information must
    be present in the matching metadata. This view omits solver factors/RHS
    and is intentionally not a GSEquilibrium suitable for a restart.
    """
    from pathlib import Path as FilePath
    from types import SimpleNamespace
    import json
    from .mesh import P1Mesh, generate_p3_nodes
    from .profiles import reconstruct_g_from_ffprime
    path = FilePath(reference_file)
    meta_path = path.parent/"metadata.json" if metadata_file is None else FilePath(metadata_file)
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    if meta.get("reference", {}).get("converged") is not True:
        raise ValueError("Matching metadata must identify a converged reference GS solve.")
    with np.load(path, allow_pickle=False) as saved:
        mesh = generate_p3_nodes(P1Mesh(saved["p1_points"], saved["p1_elements"]))
        tolerance = 64*np.finfo(float).eps*max(1., float(np.max(abs(mesh.points))))
        if (not np.array_equal(mesh.elements, saved["elements"]) or
                mesh.points.shape != saved["points"].shape or
                not np.allclose(mesh.points, saved["points"], rtol=0, atol=tolerance)):
            raise ValueError("Saved P1/P3 connectivity or node ordering is inconsistent.")
        psi = np.array(saved["psi"], dtype=float)
        axis_point = np.array(saved["axis_point"], dtype=float)
        span, axis = float(saved["dpsi"]), float(saved["psi_axis"])
        profile = reconstruct_g_from_ffprime(saved["profile_psi_n"], saved["ffprime"], float(saved["f"][-1]))
        if not np.allclose(profile.f, saved["f"], rtol=1e-12, atol=1e-12):
            raise ValueError("Saved F and FFprime profiles are inconsistent.")
    if psi.shape != (mesh.n_dof,) or not np.all(np.isfinite(psi)):
        raise ValueError("Invalid saved P3 flux vector.")
    if axis_point.shape != (2,) or not np.all(np.isfinite(axis_point)) or axis_point[0] <= 0:
        raise ValueError("Invalid saved magnetic axis.")
    if meta["reference"].get("n_dof") != mesh.n_dof or not np.isclose(
            float(meta["reference"].get("aN", np.nan)), axis_point[0], rtol=1e-12, atol=1e-12):
        raise ValueError("Reference NPZ and metadata do not match.")
    psi.setflags(write=False)
    axis_point.setflags(write=False)
    return LocalShearReference(mesh, psi, profile,
        SimpleNamespace(psi_axis=axis, dpsi=span),
        SimpleNamespace(a_minor=float(meta["definition"]["a_minor"])), axis_point, True)


def build_local_shear_map(equilibrium, *, nr=300, nz=None, hessian_method="p3",
                          poloidal_sign=1, psi_n_range=(.02, .98),
                          relative_gradient_floor=1e-8):
    """Sample a 2D local-shear map; array order is (nZ,nR), for plotting."""
    for name, value in (("nr", nr), ("nz", nz)):
        if value is not None and (isinstance(value, (bool, np.bool_)) or
                not isinstance(value, (int, np.integer)) or value < 16):
            raise ValueError(name + " must be an integer >=16, or nz=None.")
    p = equilibrium.mesh.p1_points
    r = np.linspace(p[:, 0].min(), p[:, 0].max(), nr)
    nz = nz if nz is not None else max(16, int(np.floor(nr*np.ptp(p[:, 1])/np.ptp(p[:, 0])+.5)))
    z = np.linspace(p[:, 1].min(), p[:, 1].max(), nz)
    rr, zz = np.meshgrid(r, z)
    result = compute_local_shear(equilibrium, rr, zz, hessian_method=hessian_method,
        poloidal_sign=poloidal_sign, psi_n_range=psi_n_range,
        relative_gradient_floor=relative_gradient_floor)
    if not np.any(result["valid"]):
        raise ValueError("No valid local-shear samples remain in the requested flux range.")
    arrays = (r, z, result["S"], result["psi_n"], result["valid"])
    for a in arrays:
        a.setflags(write=False)
    return LocalShearMap(*arrays, hessian_method, int(poloidal_sign),
                         tuple(psi_n_range), float(relative_gradient_floor))


def plot_local_shear(equilibrium, shear, *, levels=None, length_scale=None):
    """Return (figure, axes): negative dashed, positive solid, zero thick black.

    Raw S is used by default. length_scale=L plots L^3*S and interprets
    levels in that dimensionless normalization. No image filtering is used.
    """
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    factor = 1.
    label = r"$S$ [m$^{-3}$]"
    if length_scale is not None:
        length_scale = float(length_scale)
        if not np.isfinite(length_scale) or length_scale <= 0:
            raise ValueError("length_scale must be positive and finite.")
        factor = length_scale**3
        label = rf"$L^3 S$, $L={length_scale:g}$ m"
    values = factor*np.asarray(shear.S)
    data = np.ma.array(values, mask=~shear.valid)
    finite = values[shear.valid]
    if not finite.size or not np.all(np.isfinite(finite)):
        raise ValueError("The map has no finite valid shear data.")
    if levels is None:
        scale = float(np.percentile(abs(finite), 65))
        levels = scale*np.array([-1., -.5, -.25, 0., .25, .5, 1.])
    else:
        levels = np.asarray(levels, dtype=float)
        if levels.ndim != 1 or not len(levels) or not np.all(np.isfinite(levels)) or np.any(np.diff(levels) <= 0):
            raise ValueError("levels must be finite and strictly increasing.")
    fig, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
    lower, upper = float(finite.min()), float(finite.max())
    handles = []
    reference_scale = factor/float(equilibrium.problem.a_minor)**3
    constant = upper-lower <= 128*np.finfo(float).eps*max(reference_scale, abs(lower), abs(upper))
    if constant:
        ax.text(.5, .95, "S is numerically constant; no isolated contour",
                ha="center", va="top", transform=ax.transAxes)
    else:
        for selected, color, style, legend in (
                (levels[levels < 0], "tab:blue", "--", "S < 0"),
                (levels[levels > 0], "tab:red", "-", "S > 0")):
            selected = selected[(selected > lower) & (selected < upper)]
            if len(selected):
                cs = ax.contour(shear.r, shear.z, data, levels=selected,
                                colors=color, linestyles=style, linewidths=.85, corner_mask=False)
                ax.clabel(cs, inline=True, fontsize=7, fmt="%.2g")
                handles.append(Line2D([], [], color=color, linestyle=style, label=legend))
        if lower < 0 < upper:
            ax.contour(shear.r, shear.z, data, levels=[0.], colors="black",
                       linewidths=2.2, corner_mask=False)
            handles.append(Line2D([], [], color="black", linewidth=2.2, label="S = 0"))
        else:
            ax.text(.5, .95, "No S = 0 crossing in the plotted region", ha="center",
                    va="top", transform=ax.transAxes)
    ax.add_collection(LineCollection(equilibrium.mesh.p1_points[equilibrium.mesh.boundary_edges],
                                     colors="black", linewidths=1.1))
    ax.plot(*equilibrium.axis_point, marker="+", color="black", markersize=8, linestyle="none")
    ax.set(xlim=(shear.r[0], shear.r[-1]), ylim=(shear.z[0], shear.z[-1]),
           xlabel="R [m]", ylabel="Z [m]",
           title=f"Constant local shear ({shear.hessian_method} Hessian)\n{label}")
    ax.set_aspect("equal", adjustable="box")
    if handles:
        ax.legend(handles=handles, loc="best", fontsize=8)
    return fig, ax
