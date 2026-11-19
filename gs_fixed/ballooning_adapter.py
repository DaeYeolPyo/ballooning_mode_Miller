"""Expose a solved P3 GS equilibrium to the existing eq_process/salpha API.

Return (EquilibriumData, FEMEquilibriumFields, diagnostics), the exact input
contract accepted by salpha.prepare_fem_salpha_equilibrium(adapter=...).
The outward physical flux is Psi=abs(native_span)*psi_N. This coordinate
conversion is recorded explicitly; it is not a general COCOS conversion.
"""
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator, interp1d

from .fem import MU0
from .post import P3FieldEvaluator
from .profiles import validate_profile_grid


def _vector(values, name, size=None):
    values = np.asarray(values)
    if not np.isrealobj(values):
        raise ValueError(name + " must be real.")
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or (size is not None and values.size != size):
        raise ValueError(name + " has the wrong vector shape.")
    if not np.all(np.isfinite(values)):
        raise ValueError(name + " must be finite.")
    return values


def _grid_count(value, name):
    if (isinstance(value, (bool, np.bool_)) or
            not isinstance(value, (int, np.integer)) or value < 16):
        raise ValueError(name + " must be an integer >= 16.")
    return int(value)


def _ordered_boundary(mesh):
    """Recover one closed P1 boundary loop, with counter-clockwise order."""
    edges = np.asarray(mesh.boundary_edges)
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(int(a), []).append(int(b))
        adjacency.setdefault(int(b), []).append(int(a))
    if not adjacency or any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError("Require one closed, nonbranching mesh boundary.")
    start = min(adjacency)
    previous, current = None, start
    order, seen = [], set()
    while current not in seen:
        seen.add(current)
        order.append(current)
        neighbors = sorted(adjacency[current])
        following = neighbors[0] if neighbors[0] != previous else neighbors[1]
        previous, current = current, following
    if current != start or len(order) != len(edges):
        raise ValueError("Mesh boundary must contain exactly one closed loop.")
    points = mesh.p1_points[np.asarray(order)].copy()
    area2 = np.sum(points[:, 0]*np.roll(points[:, 1], -1)
                   - np.roll(points[:, 0], -1)*points[:, 1])
    if not np.isfinite(area2) or area2 == 0:
        raise ValueError("Degenerate mesh boundary.")
    return points if area2 > 0 else points[::-1].copy()


class FEMEquilibriumFields:
    """eq_process.EquilibriumFields-compatible P3 evaluator.

    Profile function arguments are psi_N; pprime, FFprime and qprime return
    derivatives with respect to the outward PHYSICAL flux. No derivative
    is taken from the Cartesian contour cache. Outside the FEM domain,
    spatial results are NaN. Mesh geometry is retained; no SuperLU object
    or full GS problem is retained.
    """
    def __init__(self, mesh, native_psi, native_axis, native_span, data, evaluator=None,
                 q_interpolant=None):
        self.mesh = mesh
        self._native_psi = np.array(native_psi, dtype=float, copy=True)
        self._native_psi.setflags(write=False)
        self._native_axis = float(native_axis)
        self._psi_span = abs(float(native_span))
        self._spatial_scale = self._psi_span/float(native_span)
        self._evaluator = evaluator
        for output, name in (("F", "fpol"), ("p", "pres"),
                             ("pprime", "pprime"), ("FFprime", "ffprim")):
            setattr(self, output, interp1d(
                data.psin, getattr(data, name), kind="linear", assume_sorted=True,
                bounds_error=False, fill_value=np.nan,
            ))
        self.q = (PchipInterpolator(data.psin, data.qpsi, extrapolate=True)
                  if q_interpolant is None else q_interpolant)
        self.qprimeN = self.q.derivative()

    def __getstate__(self):
        # Matplotlib's C++ triangle locator is not pickleable. Rebuild it
        # lazily in a receiving process if spatial evaluation is requested.
        state = self.__dict__.copy()
        state["_evaluator"] = None
        return state

    def _spatial(self, r, z):
        if not np.isrealobj(r) or not np.isrealobj(z):
            raise ValueError("R and Z must be real.")
        r, z = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(z, dtype=float))
        if self._evaluator is None:
            self._evaluator = P3FieldEvaluator(self.mesh)
        value, gradient, _ = self._evaluator.evaluate(
            self._native_psi, np.column_stack((r.ravel(), z.ravel())))
        value = (value-self._native_axis)*self._spatial_scale
        gradient = gradient*self._spatial_scale
        return r, z, value.reshape(r.shape), gradient.reshape(r.shape+(2,))

    def psi(self, r, z):
        return self._spatial(r, z)[2]

    def psiR(self, r, z):
        return self._spatial(r, z)[3][..., 0]

    def psiZ(self, r, z):
        return self._spatial(r, z)[3][..., 1]

    def psiN_from_psi(self, psi):
        return np.asarray(psi)/self._psi_span

    def psi_from_psiN(self, psi_n):
        return self._psi_span*np.asarray(psi_n)

    def qprime(self, psi_n):
        return self.qprimeN(psi_n)/self._psi_span

    def eval(self, r, z):
        r, z, psi, gradient = self._spatial(r, z)
        psi_n = self.psiN_from_psi(psi)
        grad_psi = np.hypot(gradient[..., 0], gradient[..., 1])
        f = self.F(psi_n)
        with np.errstate(invalid="ignore", divide="ignore"):
            bp, bphi = grad_psi/r, f/r
        return SimpleNamespace(
            R=r, Z=z, psi=psi, psiR=gradient[..., 0], psiZ=gradient[..., 1],
            psiN=psi_n, gradPsi=grad_psi, F=f, q=self.q(psi_n), p=self.p(psi_n),
            pprime=self.pprime(psi_n), ffprime=self.FFprime(psi_n),
            qprime=self.qprime(psi_n), Bp=bp, Bphi=bphi, B2=bp**2+bphi**2,
        )


def build_ballooning_adapter(equilibrium, *, nr=257, nz=None):
    """Return (data, fields, diagnostics) for the EXISTING Python Miller path.

    ``data.psirz`` has shape (nr,nz), matching eq_process. Its NaN exterior
    is intentional. Always pass the returned ``fields`` to Equilibrium,
    since EquilibriumFields(data) would rebuild grid-difference gradients.
    Profiles follow the MATLAB adapter's common-grid/PCHIP resampling.
    """
    from .solver import GSEquilibrium
    if not isinstance(equilibrium, GSEquilibrium):
        raise TypeError("equilibrium must be returned by solve_gs().")
    if not equilibrium.converged:
        raise ValueError("Refuse an unconverged GS or q-constraint solution.")
    if equilibrium.q_result is None:
        raise ValueError("The adapter requires q; solve mode 1 with compute_q=True.")
    rhs = _vector(equilibrium.state.final_rhs, "final_rhs", equilibrium.mesh.n_dof)
    return _build_adapter(equilibrium, float(np.sum(rhs)/MU0), nr=nr, nz=nz)


def load_ballooning_adapter(reference_file, *, nr=257, nz=None):
    """Adapt one converged reference_gs.npz from the C-shape scan, without GS."""
    import json
    from pathlib import Path
    from .post import load_local_shear_reference

    path = Path(reference_file)
    metadata = json.loads(path.with_name("metadata.json").read_text(encoding="utf-8"))
    reference = load_local_shear_reference(path)
    definition = metadata["definition"]
    with np.load(path, allow_pickle=False) as saved:
        equilibrium = SimpleNamespace(
            state=reference.state, mesh=reference.mesh, psi=reference.psi,
            axis_point=reference.axis_point, f_profile=reference.f_profile,
            problem=SimpleNamespace(r0=definition["r0"], b0=definition["b0"],
                                    a_minor=definition["a_minor"], psi_n=saved["profile_psi_n"]),
            pressure=saved["pressure"], pprime=saved["pprime"],
            q_result=SimpleNamespace(psi_n=saved["q_psi_n"], q=saved["q"]),
        )
    current = float(metadata["global"]["ip"])
    if not np.isfinite(current):
        raise ValueError("Invalid saved plasma current.")
    return _build_adapter(equilibrium, current, nr=nr, nz=nz, interior_q_slopes=True)


def _build_adapter(equilibrium, current, *, nr, nz, interior_q_slopes=False):
    from eq_process import EquilibriumData

    nr = _grid_count(nr, "nr")
    problem, state, mesh = equilibrium.problem, equilibrium.state, equilibrium.mesh
    psi = _vector(equilibrium.psi, "psi", mesh.n_dof)
    native_axis, native_span = float(state.psi_axis), float(state.dpsi)
    if not np.isfinite(native_axis) or not np.isfinite(native_span) or native_span == 0:
        raise ValueError("Invalid native axis or flux span.")
    span = abs(native_span)
    axis_point = np.asarray(equilibrium.axis_point, dtype=float)
    if axis_point.shape != (2,) or not np.all(np.isfinite(axis_point)) or axis_point[0] <= 0:
        raise ValueError("Invalid magnetic axis.")
    r0, b0, minor = float(problem.r0), float(problem.b0), float(problem.a_minor)
    if not np.all(np.isfinite([r0, b0, minor])) or min(r0, minor) <= 0 or b0 == 0:
        raise ValueError("Invalid R0, B0 or a_minor reference data.")
    s = validate_profile_grid(problem.psi_n)
    fp = equilibrium.f_profile
    fs = validate_profile_grid(fp.psi_n)
    if fs.shape != s.shape or not np.allclose(fs, s, rtol=0, atol=1.e-13):
        raise ValueError("F profile must use the common GS profile grid.")
    count = len(s)
    profiles = [_vector(a, name, count) for a, name in (
        (fp.f, "F"), (equilibrium.pressure, "pressure"),
        (equilibrium.pprime, "dp/dpsi_N"), (fp.ffprime, "F*dF/dpsi_N"))]
    qs = validate_profile_grid(equilibrium.q_result.psi_n)
    q = _vector(equilibrium.q_result.q, "q", len(qs))
    if np.any(q == 0) or np.any(profiles[0] == 0):
        raise ValueError("F and q must be nonzero.")
    profile_grid = np.linspace(0., 1., nr)
    fpol, pres, pprime, ffprim = [PchipInterpolator(s, a)(profile_grid) for a in profiles]
    pprime, ffprim = pprime/span, ffprim/span
    q_interpolant = None
    if interior_q_slopes:
        # q(1) is extrapolated in the GS postprocessor. Keep all saved q knots,
        # but let only the directly computed knots set the interior slopes.
        interior = PchipInterpolator(qs[:-1], q[:-1])
        edge_slope = PchipInterpolator(qs, q).derivative()(1.)
        slopes = np.r_[interior.derivative()(qs[:-1]), edge_slope]
        q_interpolant = CubicHermiteSpline(qs, q, slopes)
        qpsi = q_interpolant(profile_grid)
    else:
        # Match finalize_GS_equilibrium + build_ballooning_adapter in MATLAB.
        q_common = PchipInterpolator(qs, q)(s)
        qpsi = PchipInterpolator(s, q_common)(profile_grid)

    boundary = _ordered_boundary(mesh)
    rmin, zmin = mesh.p1_points.min(axis=0)
    rmax, zmax = mesh.p1_points.max(axis=0)
    if rmin <= 0 or rmax <= rmin or zmax <= zmin:
        raise ValueError("Invalid positive-R mesh domain.")
    nz = max(16, int(np.floor(nr*(zmax-zmin)/(rmax-rmin)+.5))) if nz is None else _grid_count(nz, "nz")
    rgrid, zgrid = np.linspace(rmin, rmax, nr), np.linspace(zmin, zmax, nz)
    rr, zz = np.meshgrid(rgrid, zgrid, indexing="ij")
    evaluator = P3FieldEvaluator(mesh)
    points = np.column_stack((rr.ravel(), zz.ravel()))
    cache = np.empty(len(points))
    # Bound temporary shape-function and sparse-operator allocations per worker.
    for start in range(0, len(points), 32768):
        block = slice(start, start + 32768)
        operator, inside = evaluator.value_operator(points[block])
        values = ((operator@psi)-native_axis)*(span/native_span)
        values[~inside] = np.nan
        cache[block] = values
    data = EquilibriumData(
        filename="", header="P3 FEM GS ballooning adapter", idum=0, nw=nr, nh=nz,
        rdim=float(rmax-rmin), zdim=float(zmax-zmin), rcentr=r0,
        rleft=float(rmin), zmid=float((zmin+zmax)/2),
        rmaxis=float(axis_point[0]), zmaxis=float(axis_point[1]),
        simag=0., sibry=span, bcentr=b0, current=current,
        fpol=fpol, pres=pres, ffprim=ffprim, pprime=pprime, qpsi=qpsi,
        psirz=cache.reshape(nr, nz), nbbbs=len(boundary), limitr=0,
        rbbbs=boundary[:, 0], zbbbs=boundary[:, 1],
        rlim=np.empty(0), zlim=np.empty(0), rgrid=rgrid, zgrid=zgrid,
        psin=profile_grid, scalar20_raw=np.empty(0),
    )
    fields = FEMEquilibriumFields(mesh, psi, native_axis, native_span, data, evaluator,
                                  q_interpolant=q_interpolant)
    diagnostics = SimpleNamespace(
        nativeAxis=native_axis, nativeBoundary=native_axis+native_span,
        nativeSpan=native_span, outwardSpan=span, spatialScale=span/native_span,
        gridSize=(nr, nz), finiteGridFraction=float(np.mean(np.isfinite(cache))),
        reference=SimpleNamespace(R0=r0, B0=b0, aMinor=minor),
        fluxConvention="Psi=abs(nativeSpan)*psi_N; native current retained as metadata",
        fieldEvaluation="direct P3 on original constrained triangles",
        qInterpolation=("PCHIP slopes from direct q knots with a Hermite edge interval"
                        if interior_q_slopes else "common-grid PCHIP"),
        profileGridSize=nr,
    )
    data.adapter = diagnostics
    return data, fields, diagnostics
