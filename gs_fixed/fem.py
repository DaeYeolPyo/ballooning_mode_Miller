"""P3 reference-element kernels for a fixed-boundary GS solver.

Geometry: affine P1 triangles.
Field: ten-node P3 Lagrange interpolation.
Dunavant weights integrate the reference triangle of area 1/2.
"""

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.interpolate import interp1d

from .mesh import P3Mesh


MU0 = 4.0 * np.pi * 1.0e-7


P3_REFERENCE_NODES = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0 / 3.0, 0.0],
    [2.0 / 3.0, 0.0],
    [2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0],
    [0.0, 2.0 / 3.0],
    [0.0, 1.0 / 3.0],
    [1.0 / 3.0, 1.0 / 3.0],
])
P3_REFERENCE_NODES.setflags(write=False)


@dataclass(frozen=True)
class P3Quadrature:
    degree: int
    xi: np.ndarray
    eta: np.ndarray
    w: np.ndarray
    N: np.ndarray
    dNdxi: np.ndarray
    dNdeta: np.ndarray

    @property
    def nq(self):
        return self.w.size


def construct_p3_shape_functions(xi, eta):
    """Return N (..., 10) and dNref (..., 2, 10).

    xi and eta may be scalars or broadcast-compatible arrays.
    Derivative rows are [d/dxi, d/deta].
    Exterior points are allowed for high-degree Dunavant rules.
    """
    xi, eta = np.broadcast_arrays(
        np.asarray(xi, dtype=float),
        np.asarray(eta, dtype=float),
    )
    if not (np.all(np.isfinite(xi)) and np.all(np.isfinite(eta))):
        raise ValueError("Reference coordinates must be finite.")

    L = np.stack((1.0 - xi - eta, xi, eta), axis=-1)
    dL = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])
    N = np.empty(xi.shape + (10,))
    dNref = np.empty(xi.shape + (2, 10))

    for i in range(3):
        x = L[..., i]
        N[..., i] = 0.5 * x * (3.0 * x - 1.0) * (3.0 * x - 2.0)
        dNref[..., :, i] = (
            (13.5 * x**2 - 9.0 * x + 1.0)[..., None] * dL[i]
        )

    # Directed edges: 1->2, 2->3, 3->1; two nodes per edge.
    pairs = ((0, 1), (1, 0), (1, 2), (2, 1), (2, 0), (0, 2))
    for k, (i, j) in enumerate(pairs, start=3):
        x, y = L[..., i], L[..., j]
        N[..., k] = 4.5 * x * y * (3.0 * x - 1.0)
        dNref[..., :, k] = (
            (4.5 * y * (6.0 * x - 1.0))[..., None] * dL[i]
            + (4.5 * x * (3.0 * x - 1.0))[..., None] * dL[j]
        )

    a, b, c = L[..., 0], L[..., 1], L[..., 2]
    N[..., 9] = 27.0 * a * b * c
    dNref[..., :, 9] = 27.0 * (
        (b * c)[..., None] * dL[0]
        + (a * c)[..., None] * dL[1]
        + (a * b)[..., None] * dL[2]
    )
    return N, dNref


def gauss_dunavant_quadrature(p=10):
    """Return xi, eta, w, each (nq,), for polynomial degree p=1..20."""
    if isinstance(p, (bool, np.bool_)) or not isinstance(p, (int, np.integer)):
        raise TypeError("p must be an integer.")
    p = int(p)
    if not 1 <= p <= 20:
        raise ValueError("Dunavant degree must be between 1 and 20.")

    w0, group1, group2 = _DUNAVANT[p]
    points, weights = [], []

    if w0 is not None:
        points.append((1.0 / 3.0, 1.0 / 3.0))
        weights.append(w0)

    for a, b, weight in group1:
        points.extend(((b, b), (a, b), (b, a)))
        weights.extend((weight,) * 3)

    for a, b, c, weight in group2:
        points.extend(((b, c), (c, b), (a, c),
                       (c, a), (a, b), (b, a)))
        weights.extend((weight,) * 6)

    points = np.asarray(points, dtype=float)
    w = 0.5 * np.asarray(weights, dtype=float)
    return points[:, 0].copy(), points[:, 1].copy(), w


def precompute_quadrature(p=10):
    """Precompute P3 basis arrays once for a selected quadrature rule."""
    xi, eta, w = gauss_dunavant_quadrature(p)
    N, dNref = construct_p3_shape_functions(xi, eta)
    dNdxi = dNref[:, 0, :].copy()
    dNdeta = dNref[:, 1, :].copy()

    for array in (xi, eta, w, N, dNdxi, dNdeta):
        array.setflags(write=False)

    return P3Quadrature(int(p), xi, eta, w, N, dNdxi, dNdeta)


def _affine_geometry(vertices):
    vertices = np.asarray(vertices, dtype=float)
    if vertices.shape != (3, 2) or not np.all(np.isfinite(vertices)):
        raise ValueError("vertices must be a finite (3, 2) array.")
    if np.any(vertices[:, 0] <= 0.0):
        raise ValueError("All vertex radii R must be positive.")

    J = np.column_stack((
        vertices[1] - vertices[0],
        vertices[2] - vertices[0],
    ))
    detJ = float(np.linalg.det(J))
    scale = np.linalg.norm(J, ord=np.inf)
    tolerance = 64.0 * np.finfo(float).eps * scale**2
    if not np.isfinite(detJ) or detJ <= tolerance:
        raise ValueError(
            "Triangle must be nondegenerate and counter-clockwise."
        )
    return vertices, J, detJ


def eval_p3_element(xi, eta, vertices):
    """Return N, dNphys, detJ, Rq, Zq on an affine physical triangle.

    N: (..., 10); dNphys: (..., 2, 10).
    Physical derivative rows are [d/dR, d/dZ].
    """
    vertices, J, detJ = _affine_geometry(vertices)
    xi, eta = np.broadcast_arrays(
        np.asarray(xi, dtype=float),
        np.asarray(eta, dtype=float),
    )
    N, dNref = construct_p3_shape_functions(xi, eta)
    dNphys = np.linalg.solve(J.T, dNref)
    Rq = vertices[0, 0] + J[0, 0] * xi + J[0, 1] * eta
    Zq = vertices[0, 1] + J[1, 0] * xi + J[1, 1] * eta

    if np.any(~np.isfinite(Rq)) or np.any(~np.isfinite(Zq)):
        raise ValueError("Physical evaluation coordinates must be finite.")
    if np.any(Rq <= 0.0):
        raise ValueError("Encountered R <= 0 at an evaluation point.")

    return N, dNphys, detJ, Rq, Zq


def construct_k_matrix(vertices, quad, sym=False):
    """Return the (10, 10) GS element stiffness matrix.

    Ke_ij = integral (1/R) grad(N_i) dot grad(N_j) dR dZ.
    Symmetry cleanup is optional and disabled by default.
    """
    if not isinstance(quad, P3Quadrature):
        raise TypeError("quad must be produced by precompute_quadrature().")

    vertices, J, detJ = _affine_geometry(vertices)
    Rq = vertices[0, 0] + J[0, 0] * quad.xi + J[0, 1] * quad.eta
    if np.any(~np.isfinite(Rq)) or np.any(Rq <= 0.0):
        raise ValueError("Quadrature-point radii must be finite and positive.")

    dNdR = (J[1, 1] * quad.dNdxi - J[1, 0] * quad.dNdeta) / detJ
    dNdZ = (-J[0, 1] * quad.dNdxi + J[0, 0] * quad.dNdeta) / detJ
    alpha = quad.w * detJ / Rq
    Ke = (
        dNdR.T @ (alpha[:, None] * dNdR)
        + dNdZ.T @ (alpha[:, None] * dNdZ)
    )
    return 0.5 * (Ke + Ke.T) if sym else Ke



def _vector(values, name, size=None):
    a = np.asarray(values)
    if np.iscomplexobj(a):
        raise ValueError(f"{name} must be real.")
    a = np.asarray(a, dtype=float)
    if a.ndim not in (1, 2) or (a.ndim == 2 and 1 not in a.shape):
        raise ValueError(f"{name} must be a vector.")
    a = a.reshape(-1)
    if size is not None and a.size != size:
        raise ValueError(f"{name} must contain {size} values.")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must contain finite values.")
    return a


def _mesh_data(mesh):
    if not isinstance(mesh, P3Mesh):
        raise TypeError("mesh must be a P3Mesh.")
    p, t, ids = mesh.p1_points, mesh.p1_elements, mesh.elements
    ndof = mesh.n_dof
    if (
        mesh.points.shape != (ndof, 2)
        or p.ndim != 2 or p.shape[1] != 2
        or ids.ndim != 2 or ids.shape[1] != 10 or len(ids) < 1
        or t.shape != (len(ids), 3)
    ):
        raise ValueError("Incompatible P3 mesh array shapes.")
    if not (
        np.all(np.isfinite(mesh.points)) and np.all(np.isfinite(p))
        and np.all(mesh.points[:, 0] > 0) and np.all(p[:, 0] > 0)
    ):
        raise ValueError("Mesh coordinates must be finite with R > 0.")
    for values, bound in ((t, len(p)), (ids, ndof)):
        if (
            not np.issubdtype(values.dtype, np.integer)
            or np.any(values < 0) or np.any(values >= bound)
        ):
            raise ValueError("Invalid mesh connectivity.")
    if np.any(np.diff(np.sort(ids, axis=1), axis=1) == 0):
        raise ValueError("An element contains duplicate P3 node IDs.")
    if not np.array_equal(ids[:, :3], t):
        raise ValueError("P1/P3 vertex connectivity must agree.")
    return p, t, ids, ndof


def assemble_k_matrix(mesh, quad, sym=False):
    """Assemble the global GS stiffness matrix as scipy CSR."""
    p, t, connectivity, ndof = _mesh_data(mesh)
    if not isinstance(quad, P3Quadrature):
        raise TypeError("quad must be a P3Quadrature.")

    nentries = 100 * len(connectivity)
    rows = np.empty(nentries, dtype=np.int64)
    cols = np.empty(nentries, dtype=np.int64)
    values = np.empty(nentries)

    for e, ids in enumerate(connectivity):
        ke = construct_k_matrix(p[t[e]], quad, sym=sym)
        block = slice(100 * e, 100 * (e + 1))
        rows[block] = np.repeat(ids, 10)
        cols[block] = np.tile(ids, 10)
        values[block] = ke.ravel()

    if not np.all(np.isfinite(values)):
        raise ValueError("The stiffness matrix contains NaN or Inf.")
    return sparse.coo_matrix(
        (values, (rows, cols)), shape=(ndof, ndof)
    ).tocsr()


def _prepare_rhs_profiles(pprime, ffprime, dpsi):
    pprime = _vector(pprime, "pprime")
    ffprime = _vector(ffprime, "ffprime", pprime.size)
    if pprime.size < 2:
        raise ValueError("At least two profile samples are required.")

    span = np.asarray(dpsi)
    if (
        span.ndim != 0 or np.iscomplexobj(span)
        or np.issubdtype(span.dtype, np.bool_)
    ):
        raise ValueError("dpsi must be a real nonzero scalar.")
    span = float(span)
    if not np.isfinite(span) or span == 0:
        raise ValueError("dpsi must be a finite nonzero scalar.")

    grid = np.linspace(0.0, 1.0, pprime.size)
    physical_profiles = np.column_stack((pprime, ffprime)) / span
    if not np.all(np.isfinite(physical_profiles)):
        raise ValueError("Profile conversion produced NaN or Inf.")

    profiles = interp1d(
        grid, physical_profiles, axis=0, kind="linear",
        bounds_error=False, fill_value="extrapolate",
        assume_sorted=True,
    )
    return span, profiles


def _construct_rhs_prepared(vertices, psi_e, quad, dpsi, profiles):
    vertices, J, detJ = _affine_geometry(vertices)
    rq = vertices[0, 0] + J[0, 0] * quad.xi + J[0, 1] * quad.eta
    if np.any(~np.isfinite(rq)) or np.any(rq <= 0):
        raise ValueError("Quadrature radii must be finite and positive.")

    psi_q = quad.N @ psi_e
    psi_n = (psi_q + dpsi) / dpsi
    if not np.all(np.isfinite(psi_n)):
        raise ValueError("Normalized flux contains NaN or Inf.")

    physical = profiles(psi_n)
    source = MU0 * rq * physical[:, 0] + physical[:, 1] / rq
    fe = quad.N.T @ (quad.w * detJ * source)
    if not np.all(np.isfinite(fe)):
        raise ValueError("The element RHS contains NaN or Inf.")
    return fe


def construct_rhs(vertices, psi_e, quad, pprime, ffprime, dpsi):
    """Return a (10,) element GS load vector.

    pprime:  dp/dpsi_N on a uniform grid spanning [0, 1].
    ffprime: F*dF/dpsi_N on the same grid.
    dpsi:    signed psi_boundary - psi_axis.

    This RHS follows the current MATLAB convention psi_boundary = 0.
    Thus psi_N = (psi + dpsi)/dpsi.

    Linear interpolation AND linear extrapolation are used.
    Flux values are not clipped to [0, 1].
    """
    if not isinstance(quad, P3Quadrature):
        raise TypeError("quad must be a P3Quadrature.")
    psi_e = _vector(psi_e, "psi_e", 10)
    span, profiles = _prepare_rhs_profiles(pprime, ffprime, dpsi)
    return _construct_rhs_prepared(vertices, psi_e, quad, span, profiles)


def assemble_rhs(mesh, psi, quad, pprime, ffprime, dpsi):
    """Assemble the global (ndof,) nonlinear GS RHS.

    Profile interpolants are constructed once per global assembly.
    Requires the same zero-boundary-flux convention as construct_rhs().
    """
    p, t, connectivity, ndof = _mesh_data(mesh)
    if not isinstance(quad, P3Quadrature):
        raise TypeError("quad must be a P3Quadrature.")
    psi = _vector(psi, "psi", ndof)
    span, profiles = _prepare_rhs_profiles(pprime, ffprime, dpsi)

    rhs = np.zeros(ndof)
    for e, ids in enumerate(connectivity):
        fe = _construct_rhs_prepared(
            p[t[e]], psi[ids], quad, span, profiles
        )
        rhs[ids] += fe
    if not np.all(np.isfinite(rhs)):
        raise ValueError("The assembled RHS contains NaN or Inf.")
    return rhs


@dataclass(frozen=True)
class DirichletBC:
    """Reusable partition and RHS correction for fixed Dirichlet data."""
    n_dof: int
    free_nodes: np.ndarray
    boundary_nodes: np.ndarray
    boundary_values: np.ndarray
    boundary_correction: np.ndarray

    def reduce_rhs(self, rhs):
        """Apply the cached correction to a new global RHS."""
        rhs = _vector(rhs, "rhs", self.n_dof)
        return rhs[self.free_nodes] - self.boundary_correction

    def restore(self, psi_free):
        """Restore a full nodal vector from the reduced solution."""
        psi_free = _vector(psi_free, "psi_free", self.free_nodes.size)
        psi = np.empty(self.n_dof)
        psi[self.free_nodes] = psi_free
        psi[self.boundary_nodes] = self.boundary_values
        return psi


def apply_bc(mesh, K, rhs, psi_boundary=0.0):
    """Return KII, rhsI, bc for KII*psiI = rhsI.

    psi_boundary accepts:
      - a scalar;
      - a vector ordered according to mesh.boundary_nodes;
      - a full nodal vector;
      - a callable psi_boundary(R, Z) evaluated at boundary nodes.

    The linear boundary reduction supports nonzero Dirichlet data.
    The GS RHS above retains its psi_boundary = 0 convention.
    """
    _, _, _, ndof = _mesh_data(mesh)
    if not sparse.issparse(K) or K.shape != (ndof, ndof):
        raise ValueError("K must be a sparse (ndof, ndof) matrix.")
    K = K.tocsr()
    if np.iscomplexobj(K.data) or not np.all(np.isfinite(K.data)):
        raise ValueError("K must contain finite real values.")
    K = K.astype(float, copy=False).tocsr()
    rhs = _vector(rhs, "rhs", ndof)

    boundary = np.asarray(mesh.boundary_nodes)
    if (
        boundary.ndim != 1 or not boundary.size
        or not np.issubdtype(boundary.dtype, np.integer)
        or np.any(boundary < 0) or np.any(boundary >= ndof)
        or np.unique(boundary).size != boundary.size
    ):
        raise ValueError("Invalid boundary-node indices.")
    boundary = boundary.copy()
    mask = np.ones(ndof, dtype=bool)
    mask[boundary] = False
    free = np.flatnonzero(mask)
    if not free.size:
        raise ValueError("No free DOFs remain.")

    if callable(psi_boundary):
        x = mesh.points[boundary]
        prescribed = psi_boundary(x[:, 0], x[:, 1])
    else:
        prescribed = psi_boundary
    prescribed = np.asarray(prescribed)
    if np.iscomplexobj(prescribed):
        raise ValueError("Boundary values must be real.")

    if prescribed.ndim == 0:
        value = float(prescribed)
        if not np.isfinite(value):
            raise ValueError("Boundary values must be finite.")
        values = np.full(boundary.size, value)
    else:
        prescribed = _vector(prescribed, "psi_boundary")
        if prescribed.size == boundary.size:
            values = prescribed.copy()
        elif prescribed.size == ndof:
            values = prescribed[boundary].copy()
        else:
            raise ValueError("Boundary-value vector has incompatible length.")

    KII = K[free, :][:, free].tocsr()
    correction = np.asarray(K[free, :][:, boundary] @ values).reshape(-1)
    if not np.all(np.isfinite(correction)):
        raise ValueError("Boundary correction contains NaN or Inf.")

    for array in (free, boundary, values, correction):
        array.setflags(write=False)
    bc = DirichletBC(ndof, free, boundary, values, correction)
    return KII, bc.reduce_rhs(rhs), bc


# Dunavant orbit tables; stored weights sum to 1.
_DUNAVANT = {
    1: (1.0, [
    ], [
    ]),
    2: (None, [
        (0.666666666666667, 0.166666666666667, 0.333333333333333),
    ], [
    ]),
    3: (-0.5625, [
        (0.6, 0.2, 0.520833333333333),
    ], [
    ]),
    4: (None, [
        (0.10810301816807, 0.445948490915965, 0.223381589678011),
        (0.816847572980459, 0.091576213509771, 0.109951743655322),
    ], [
    ]),
    5: (0.225, [
        (0.05971587178977, 0.470142064105115, 0.132394152788506),
        (0.797426985353087, 0.101286507323456, 0.125939180544827),
    ], [
    ]),
    6: (None, [
        (0.501426509658179, 0.24928674517091, 0.116786275726379),
        (0.873821971016996, 0.063089014491502, 0.050844906370207),
    ], [
        (0.053145049844817, 0.310352451033784, 0.636502499121399, 0.082851075618374),
    ]),
    7: (-0.149570044467682, [
        (0.47930806784192, 0.26034596607904, 0.175615257433208),
        (0.869739794195568, 0.065130102902216, 0.053347235608838),
    ], [
        (0.048690315425316, 0.312865496004874, 0.63844418856981, 0.077113760890257),
    ]),
    8: (0.144315607677787, [
        (0.081414823414554, 0.459292588292723, 0.095091634267285),
        (0.65886138449648, 0.17056930775176, 0.103217370534718),
        (0.898905543365938, 0.050547228317031, 0.032458497623198),
    ], [
        (0.008394777409958, 0.263112829634638, 0.728492392955404, 0.027230314174435),
    ]),
    9: (0.097135796282799, [
        (0.020634961602525, 0.489682519198738, 0.031334700227139),
        (0.125820817014127, 0.437089591492937, 0.077827541004774),
        (0.623592928761935, 0.188203535619033, 0.07964773892721),
        (0.910540973211095, 0.044729513394453, 0.025577675658698),
    ], [
        (0.036838412054736, 0.221962989160766, 0.741198598784498, 0.043283539377289),
    ]),
    10: (0.090817990382754, [
        (0.028844733232685, 0.485577633383657, 0.036725957756467),
        (0.781036849029926, 0.109481575485037, 0.045321059435528),
    ], [
        (0.14170721941488, 0.307939838764121, 0.550352941820999, 0.07275791684542),
        (0.025003534762686, 0.246672560639903, 0.728323904597411, 0.028327242531057),
        (0.009540815400299, 0.0668032510122, 0.9236559335875, 0.009421666963733),
    ]),
    11: (None, [
        (-0.069222096541517, 0.534611048270758, 0.000927006328961),
        (0.20206139406829, 0.398969302965855, 0.077149534914813),
        (0.593380199137435, 0.203309900431282, 0.059322977380774),
        (0.761298175434837, 0.119350912282581, 0.036184540503418),
        (0.935270103777448, 0.032364948111276, 0.013659731002678),
    ], [
        (0.050178138310495, 0.356620648261293, 0.593201213428213, 0.052337111962204),
        (0.021022016536166, 0.171488980304042, 0.807489003159792, 0.020707659639141),
    ]),
    12: (None, [
        (0.02356522045239, 0.488217389773805, 0.025731066440455),
        (0.120551215411079, 0.43972439229446, 0.043692544538038),
        (0.457579229975768, 0.271210385012116, 0.062858224217885),
        (0.744847708916828, 0.127576145541586, 0.034796112930709),
        (0.957365299093579, 0.02131735045321, 0.006166261051559),
    ], [
        (0.115343494534698, 0.275713269685514, 0.608943235779788, 0.040371557766381),
        (0.022838332222257, 0.28132558098994, 0.695836086787803, 0.022356773202303),
        (0.02573405054833, 0.116251915907597, 0.858014033544073, 0.017316231108659),
    ]),
    13: (0.052520923400802, [
        (0.009903630120591, 0.495048184939705, 0.01128014520933),
        (0.062566729780852, 0.468716635109574, 0.031423518362454),
        (0.170957326397447, 0.414521336801277, 0.047072502504194),
        (0.541200855914337, 0.229399572042831, 0.047363586536355),
        (0.77115100960734, 0.11442449519633, 0.031167529045794),
        (0.950377217273082, 0.024811391363459, 0.007975771465074),
    ], [
        (0.094853828379579, 0.268794997058761, 0.63635117456166, 0.036848402728732),
        (0.018100773278807, 0.291730066734288, 0.690169159986905, 0.017401463303822),
        (0.02223307667409, 0.126357385491669, 0.851409537834241, 0.015521786839045),
    ]),
    14: (None, [
        (0.022072179275643, 0.488963910362179, 0.021883581369429),
        (0.164710561319092, 0.417644719340454, 0.032788353544125),
        (0.453044943382323, 0.273477528308839, 0.051774104507292),
        (0.645588935174913, 0.177205532412543, 0.042162588736993),
        (0.876400233818255, 0.061799883090873, 0.014433699669777),
        (0.961218077502598, 0.019390961248701, 0.0049234036024),
    ], [
        (0.057124757403648, 0.172266687821356, 0.770608554774996, 0.024665753212564),
        (0.092916249356972, 0.336861459796345, 0.570222290846683, 0.038571510787061),
        (0.014646950055654, 0.298372882136258, 0.686980167808088, 0.014436308113534),
        (0.001268330932872, 0.118974497696957, 0.879757171370171, 0.005010228838501),
    ]),
    15: (None, [
        (-0.013945833716486, 0.506972916858243, 0.001916875642849),
        (0.137187291433955, 0.431406354283023, 0.044249027271145),
        (0.444612710305711, 0.277693644847144, 0.051186548718852),
        (0.747070217917492, 0.126464891041254, 0.023687735870688),
        (0.858383228050628, 0.070808385974686, 0.013289775690021),
        (0.962069659517853, 0.018965170241073, 0.004748916608192),
    ], [
        (0.133734161966621, 0.261311371140087, 0.604954466893291, 0.038550072599593),
        (0.036366677396917, 0.388046767090269, 0.575586555512814, 0.027215814320624),
        (-0.010174883126571, 0.285712220049916, 0.724462663076655, 0.002182077366797),
        (0.036843869875878, 0.215599664072284, 0.747556466051838, 0.021505319847731),
        (0.012459809331199, 0.103575616576386, 0.883964574092416, 0.007673942631049),
    ]),
    16: (0.046875697427642, [
        (0.005238916103123, 0.497380541948438, 0.006405878578585),
        (0.173061122901295, 0.413469438549352, 0.041710296739387),
        (0.059082801866017, 0.470458599066991, 0.026891484250064),
        (0.518892500060958, 0.240553749969521, 0.04213252276165),
        (0.704068411554854, 0.147965794222573, 0.030000266842773),
        (0.849069624685052, 0.075465187657474, 0.014200098925024),
        (0.96680719475395, 0.016596402623025, 0.003582462351273),
    ], [
        (0.103575692245252, 0.296555596579887, 0.599868711174861, 0.032773147460627),
        (0.020083411655416, 0.337723063403079, 0.642193524941505, 0.015298306248441),
        (-0.004341002614139, 0.204748281642812, 0.799592720971327, 0.002386244192839),
        (0.04194178646801, 0.189358492130623, 0.768699721401368, 0.019084792755899),
        (0.014317320230681, 0.085283615682657, 0.900399064086661, 0.006850054546542),
    ]),
    17: (0.033437199290803, [
        (0.005658918886452, 0.497170540556774, 0.005093415440507),
        (0.035647354750751, 0.482176322624625, 0.014670864527638),
        (0.099520061958437, 0.450239969020782, 0.024350878353672),
        (0.199467521245206, 0.400266239377397, 0.031107550868969),
        (0.495717464058095, 0.252141267970953, 0.03125711121862),
        (0.675905990683077, 0.162047004658461, 0.024815654339665),
        (0.848248235478508, 0.075875882260746, 0.014056073070557),
        (0.968690546064356, 0.015654726967822, 0.003194676173779),
    ], [
        (0.010186928826919, 0.334319867363658, 0.655493203809423, 0.008119655318993),
        (0.135440871671036, 0.292221537796944, 0.57233759053202, 0.026805742283163),
        (0.054423924290583, 0.31957488542319, 0.626001190286228, 0.018459993210822),
        (0.012868560833637, 0.190704224192292, 0.796427214974071, 0.008476868534328),
        (0.067165782413524, 0.180483211648746, 0.752351005937729, 0.018292796770025),
        (0.014663182224828, 0.080711313679564, 0.904625504095608, 0.006665632004165),
    ]),
    18: (0.030809939937647, [
        (0.013310382738157, 0.493344808630921, 0.009072436679404),
        (0.061578811516086, 0.469210594241957, 0.018761316939594),
        (0.127437208225989, 0.436281395887006, 0.019441097985477),
        (0.210307658653168, 0.394846170673416, 0.02775394861081),
        (0.500410862393686, 0.249794568803157, 0.032256225351457),
        (0.677135612512315, 0.161432193743843, 0.025074032616922),
        (0.846803545029257, 0.076598227485371, 0.015271927971832),
        (0.9514951212931, 0.02425243935345, 0.006793922022963),
        (0.913707265566071, 0.043146367216965, -0.00222309872992),
    ], [
        (0.00843053620242, 0.358911494940944, 0.632657968856636, 0.006331914076406),
        (0.131186551737188, 0.294402476751957, 0.574410971510855, 0.027257538049138),
        (0.050203151565675, 0.325017801641814, 0.624779046792512, 0.017676785649465),
        (0.066329263810916, 0.184737559666046, 0.748933176523037, 0.01837948463807),
        (0.011996194566236, 0.218796800013321, 0.769207005420443, 0.008104732808192),
        (0.014858100590125, 0.101179597136408, 0.883962302273467, 0.007634129070725),
        (-0.035222015287949, 0.020874755282586, 1.014347260005363, 4.6187660794e-05),
    ]),
    19: (0.032906331388919, [
        (0.020780025853987, 0.489609987073006, 0.010330731891272),
        (0.090926214604215, 0.454536892697893, 0.022387247263016),
        (0.197166638701138, 0.401416680649431, 0.030266125869468),
        (0.488896691193805, 0.255551654403098, 0.030490967802198),
        (0.645844115695741, 0.17707794215213, 0.024159212741641),
        (0.779877893544096, 0.110061053227952, 0.016050803586801),
        (0.888942751496321, 0.05552862425184, 0.008084580261784),
        (0.974756272445543, 0.012621863777229, 0.002079362027485),
    ], [
        (0.003611417848412, 0.395754787356943, 0.600633794794645, 0.003884876904981),
        (0.13446675453078, 0.307929983880436, 0.557603261588784, 0.025574160612022),
        (0.014446025776115, 0.26456694840652, 0.720987025817365, 0.008880903573338),
        (0.046933578838178, 0.358539352205951, 0.594527068955871, 0.016124546761731),
        (0.002861120350567, 0.157807405968595, 0.839331473680839, 0.002491941817491),
        (0.223861424097916, 0.075050596975911, 0.701087978926173, 0.018242840118951),
        (0.03464707481676, 0.142421601113383, 0.822931324069857, 0.010258563736199),
        (0.010161119296278, 0.065494628082938, 0.924344252620784, 0.003799928855302),
    ]),
    20: (0.033057055541624, [
        (-0.0019009287044, 0.5009504643522, 0.000867019185663),
        (0.023574084130543, 0.488212957934729, 0.011660052716448),
        (0.089726636099435, 0.455136681950283, 0.022876936356421),
        (0.196007481363421, 0.401996259318289, 0.030448982673938),
        (0.488214180481157, 0.255892909759421, 0.030624891725355),
        (0.647023488009788, 0.176488255995106, 0.0243680576768),
        (0.791658289326483, 0.104170855336758, 0.015997432032024),
        (0.89386207231814, 0.05306896384093, 0.007698301815602),
        (0.916762569607942, 0.041618715196029, -0.000632060497488),
        (0.976836157186356, 0.011581921406822, 0.001751134301193),
    ], [
        (0.048741583664839, 0.344855770229001, 0.60640264610616, 0.016465839189576),
        (0.006314115948605, 0.377843269594854, 0.615842614456541, 0.004839033540485),
        (0.134316520547348, 0.306635479062357, 0.559048000390295, 0.02580490653465),
        (0.013973893962392, 0.249419362774742, 0.736606743262866, 0.008471091054441),
        (0.075549132909764, 0.212775724802802, 0.711675142287434, 0.01835491410628),
        (-0.008368153208227, 0.146965436053239, 0.861402717154987, 0.000704404677908),
        (0.026686063258714, 0.137726978828923, 0.835586957912363, 0.010112684927462),
        (0.010547719294141, 0.059696109149007, 0.929756171556853, 0.00357390938595),
    ]),
}
