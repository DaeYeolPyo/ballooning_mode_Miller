"""Affine P1 meshes and ten-node P3 connectivity for fixed-boundary GS.

All connectivity uses zero-based indices.
Boundary edges are canonical pairs (smaller vertex ID first).
The meshing interface supports one simple closed boundary, without holes.
"""

from dataclasses import dataclass, field

import numpy as np


def _readonly(array):
    array.setflags(write=False)
    return array


def _points(values, name="points"):
    p = np.array(values, dtype=float, copy=True)
    if p.ndim != 2 or p.shape[1] != 2 or len(p) < 1:
        raise ValueError(f"{name} must have shape (n, 2), n >= 1.")
    if not np.all(np.isfinite(p)) or np.any(p[:, 0] <= 0):
        raise ValueError(f"{name} must be finite with R > 0.")
    return p


def _indices(values, width, npoints, name):
    t = np.asarray(values)
    if (
        t.ndim != 2 or t.shape[1] != width or len(t) < 1
        or not np.issubdtype(t.dtype, np.integer)
    ):
        raise ValueError(f"{name} must be a nonempty integer (n, {width}) array.")
    if np.any(t < 0) or np.any(t >= npoints):
        raise ValueError(f"{name} contains out-of-range zero-based indices.")
    return np.array(t, dtype=np.int64, copy=True)


def _cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def _coordinate_tolerance(p):
    return 128.0 * np.finfo(float).eps * np.max(np.abs(p))


def _edge_topology(t):
    """Stable edge order matches MATLAB: all 12, then 23, then 31."""
    nt = len(t)
    directed = np.concatenate((t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]))
    canonical = np.sort(directed, axis=1)
    edges, first, inverse, counts = np.unique(
        canonical, axis=0, return_index=True,
        return_inverse=True, return_counts=True,
    )
    order = np.argsort(first)
    remap = np.empty(len(order), dtype=np.int64)
    remap[order] = np.arange(len(order))
    ids = remap[inverse]
    edges, counts = edges[order], counts[order]

    if np.any(counts > 2):
        raise ValueError("A mesh edge belongs to more than two triangles.")
    signs = np.where(directed[:, 0] < directed[:, 1], 1.0, -1.0)
    signed_counts = np.bincount(ids, weights=signs, minlength=len(edges))
    if np.any((counts == 2) & (signed_counts != 0)):
        raise ValueError("Adjacent triangles traverse a shared edge equally.")

    element_edges = ids.reshape(3, nt).T.copy()
    boundary_edge_ids = np.flatnonzero(counts == 1)
    boundary_edges = edges[boundary_edge_ids].copy()
    if len(boundary_edges) == 0:
        raise ValueError("Mesh has no boundary edges.")
    degrees = np.bincount(boundary_edges.ravel(), minlength=t.max() + 1)
    if np.any((degrees != 0) & (degrees != 2)):
        raise ValueError("Boundary edges must form disjoint simple loops.")
    return edges, element_edges, boundary_edge_ids, boundary_edges


@dataclass(frozen=True)
class P1Mesh:
    points: np.ndarray
    elements: np.ndarray
    edges: np.ndarray = field(init=False, repr=False)
    element_edges: np.ndarray = field(init=False, repr=False)
    boundary_edge_ids: np.ndarray = field(init=False, repr=False)
    boundary_edges: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        p = _points(self.points)
        t = _indices(self.elements, 3, len(p), "elements")
        if len(np.unique(p, axis=0)) != len(p):
            raise ValueError("Mesh contains duplicate point coordinates.")
        if len(np.unique(t)) != len(p):
            raise ValueError("Mesh contains unused points.")
        if len(np.unique(np.sort(t, axis=1), axis=0)) != len(t):
            raise ValueError("Mesh contains duplicate triangles.")

        e1 = p[t[:, 1]] - p[t[:, 0]]
        e2 = p[t[:, 2]] - p[t[:, 0]]
        det = _cross(e1, e2)
        scale = np.max(np.abs(e1) + np.abs(e2), axis=1)
        tol = 64.0 * np.finfo(float).eps * scale**2
        if np.any(~np.isfinite(det)) or np.any(np.abs(det) <= tol):
            raise ValueError("Mesh contains a degenerate triangle.")
        clockwise = det < 0
        t[clockwise] = t[clockwise][:, [0, 2, 1]]

        arrays = dict(
            points=p, elements=t,
            **dict(zip(
                ("edges", "element_edges", "boundary_edge_ids", "boundary_edges"),
                _edge_topology(t),
            )),
        )
        for name, value in arrays.items():
            object.__setattr__(self, name, _readonly(value))

    @property
    def n_vertices(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.elements)


@dataclass(frozen=True)
class P3Mesh:
    points: np.ndarray
    elements: np.ndarray
    vertices: np.ndarray
    edges: np.ndarray
    edge_nodes: np.ndarray
    interior_nodes: np.ndarray
    element_edges: np.ndarray
    boundary_edges: np.ndarray
    boundary_nodes: np.ndarray
    p1_points: np.ndarray
    p1_elements: np.ndarray

    @property
    def n_dof(self):
        return len(self.points)

    @property
    def n_elements(self):
        return len(self.elements)

    @property
    def n_vertices(self):
        return len(self.p1_points)

    @property
    def n_edges(self):
        return len(self.edges)

    @property
    def free_nodes(self):
        return np.setdiff1d(
            np.arange(self.n_dof), self.boundary_nodes, assume_unique=True
        )


def _boundary_loop(points, segments):
    neighbors = {}
    for a, b in segments:
        if a == b:
            raise ValueError("A boundary segment has identical endpoints.")
        neighbors.setdefault(int(a), []).append(int(b))
        neighbors.setdefault(int(b), []).append(int(a))
    if any(len(v) != 2 or v[0] == v[1] for v in neighbors.values()):
        raise ValueError("Boundary segments must form one simple closed loop.")
    start = int(segments[0, 0])
    order, seen = [], set()
    previous, current = None, start
    while current not in seen:
        seen.add(current)
        order.append(current)
        a, b = neighbors[current]
        previous, current = current, b if a == previous else a
    if current != start or len(order) != len(neighbors):
        raise ValueError("Boundary must contain exactly one closed loop.")
    return points[order]


def _polygon_area(p):
    q = p - p[0]
    return 0.5 * abs(np.sum(_cross(q, np.roll(q, -1, axis=0))))


def _validate_polygon(p):
    """Reject crossing, touching and overlapping nonadjacent segments."""
    if len(p) < 3:
        raise ValueError("At least three boundary vertices are required.")
    tol = _coordinate_tolerance(p)
    edges = np.roll(p, -1, axis=0) - p
    lengths = np.linalg.norm(edges, axis=1)
    if np.any(lengths <= tol):
        raise ValueError("Boundary has a zero-length or unresolved segment.")

    previous = np.roll(edges, 1, axis=0)
    collinear = np.abs(_cross(previous, edges)) <= (
        tol * (lengths + np.roll(lengths, 1))
    )
    if np.any(collinear & (np.sum(previous * edges, axis=1) < 0)):
        raise ValueError("Boundary contains adjacent overlapping segments.")

    n = len(p)
    for i in range(n):
        j = np.arange(i + 2, n)
        if i == 0:
            j = j[j != n - 1]
        if not len(j):
            continue
        a, b = p[i], p[(i + 1) % n]
        c, d = p[j], p[(j + 1) % n]
        overlap = np.all(
            np.maximum(np.minimum(a, b), np.minimum(c, d))
            <= np.minimum(np.maximum(a, b), np.maximum(c, d)) + tol,
            axis=1,
        )
        ab, cd = b - a, d - c
        o1, o2 = _cross(ab, c - a), _cross(ab, d - a)
        o3, o4 = _cross(cd, a - c), _cross(cd, b - c)
        t1, t2 = tol * lengths[i], tol * lengths[j]
        intersects = (
            overlap
            & (np.minimum(o1, o2) <= t1)
            & (np.maximum(o1, o2) >= -t1)
            & (np.minimum(o3, o4) <= t2)
            & (np.maximum(o3, o4) >= -t2)
        )
        if np.any(intersects):
            raise ValueError("Boundary must not intersect or touch itself.")

    area = _polygon_area(p)
    if area <= tol * np.sum(lengths):
        raise ValueError("Boundary encloses negligible area.")


def _classify_points(points, polygon):
    """Return ray-casting inside mask and a separate boundary mask."""
    x, y = points[:, 0], points[:, 1]
    inside = np.zeros(len(points), dtype=bool)
    on_boundary = np.zeros(len(points), dtype=bool)
    tol = _coordinate_tolerance(polygon)

    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        edge = b - a
        rel = points - a
        length = np.linalg.norm(edge)
        projection = rel @ edge
        on_boundary |= (
            (np.abs(_cross(edge, rel)) <= tol * length)
            & (projection >= -tol * length)
            & (projection <= length**2 + tol * length)
        )
        if edge[1] != 0:
            crosses = (a[1] > y) != (b[1] > y)
            intersection_x = a[0] + (y - a[1]) * edge[0] / edge[1]
            inside ^= crosses & (x < intersection_x)
    return inside, on_boundary


def generate_p1_mesh(vertices, boundary, *, max_area=None):
    """Constrained triangulation of one simple closed polygon.

    vertices: (nv, 2) coordinates [R, Z], including optional interior points.
    boundary: (nb, 2) zero-based segment endpoint IDs.
    max_area: optional positive triangle-area bound [m^2].

    Triangle may subdivide segments when refinement is requested.
    Returned topology always uses the actual triangulation.
    """
    p = _points(vertices, "vertices")
    s = _indices(boundary, 2, len(p), "boundary")
    if len(np.unique(p, axis=0)) != len(p):
        raise ValueError("vertices contains duplicate coordinates.")
    polygon = _boundary_loop(p, s)
    _validate_polygon(polygon)
    inside, on_boundary = _classify_points(p, polygon)
    if not np.all(inside | on_boundary):
        raise ValueError("All supplied points must lie inside the boundary.")

    options = "pQ"
    if max_area is not None:
        max_area = float(max_area)
        if not np.isfinite(max_area) or max_area <= 0:
            raise ValueError("max_area must be finite and positive.")
        # Triangle parses decimal option values; avoid scientific notation.
        options += "a" + np.format_float_positional(max_area, trim="-")

    try:
        import triangle
    except ImportError as exc:
        raise ImportError(
            "Constrained meshing requires: python -m pip install triangle"
        ) from exc

    result = triangle.triangulate(
        {"vertices": p.copy(), "segments": s.copy()}, options
    )
    if "triangles" not in result or "segments" not in result:
        raise RuntimeError("Triangle failed to generate a constrained mesh.")
    mesh = P1Mesh(result["vertices"], result["triangles"])

    returned_segments = {
        tuple(sorted(edge)) for edge in result["segments"]
    }
    if returned_segments != {tuple(edge) for edge in mesh.boundary_edges}:
        raise RuntimeError("Triangulation boundary does not match its constraints.")

    centers = mesh.points[mesh.elements].mean(axis=1)
    inside, on_boundary = _classify_points(centers, polygon)
    if not np.all(inside | on_boundary):
        raise RuntimeError("Triangulation contains an exterior triangle.")
    v = mesh.points[mesh.elements]
    area = 0.5 * np.sum(_cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]))
    if not np.isclose(area, _polygon_area(polygon), rtol=1e-10, atol=0):
        raise RuntimeError("Triangulation does not cover the prescribed domain.")
    return mesh


def build_p1_mesh(boundary_points, *, nr=55, nz=65,
                  extra_points=None, max_area=None):
    """Build a P1 mesh from an ordered [R, Z] boundary and an interior grid.

    A repeated closing endpoint is accepted and removed.
    nr/nz specify candidate Cartesian grid sizes, not mesh node counts.
    """
    polygon = _points(boundary_points, "boundary_points")
    tol = _coordinate_tolerance(polygon)
    if len(polygon) > 1 and np.linalg.norm(polygon[-1] - polygon[0]) <= tol:
        polygon = polygon[:-1].copy()
    _validate_polygon(polygon)

    for name, value in (("nr", nr), ("nz", nz)):
        if (
            isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value < 2
        ):
            raise ValueError(f"{name} must be an integer >= 2.")

    r = np.linspace(polygon[:, 0].min(), polygon[:, 0].max(), nr)
    z = np.linspace(polygon[:, 1].min(), polygon[:, 1].max(), nz)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    grid = np.column_stack((rr.ravel(), zz.ravel()))
    inside, on_boundary = _classify_points(grid, polygon)
    grid = grid[inside & ~on_boundary]

    extra = np.empty((0, 2))
    if extra_points is not None:
        extra = _points(extra_points, "extra_points")
        # Retain supplied coordinates when an extra point matches the grid.
        for point in extra:
            grid = grid[np.linalg.norm(grid - point, axis=1) > tol]

    points = np.vstack((polygon, grid, extra))
    _, first = np.unique(points, axis=0, return_index=True)
    points = points[np.sort(first)]
    n = len(polygon)
    segments = np.column_stack((np.arange(n), np.roll(np.arange(n), -1)))
    return generate_p1_mesh(points, segments, max_area=max_area)


def generate_p3_nodes(mesh_p1):
    """Upgrade a P1Mesh to affine-geometry, cubic-field P3 connectivity."""
    if not isinstance(mesh_p1, P1Mesh):
        raise TypeError("mesh_p1 must be a P1Mesh.")

    p, t = mesh_p1.points, mesh_p1.elements
    edges = mesh_p1.edges
    nv, nt, ne = len(p), len(t), len(edges)
    edge_nodes = np.column_stack((
        nv + np.arange(ne), nv + ne + np.arange(ne),
    ))
    interior_nodes = nv + 2 * ne + np.arange(nt)

    a, b = p[edges[:, 0]], p[edges[:, 1]]
    points = np.vstack((
        p, (2 * a + b) / 3, (a + 2 * b) / 3, p[t].mean(axis=1),
    ))
    elements = np.empty((nt, 10), dtype=np.int64)
    elements[:, :3] = t
    for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
        pair = edge_nodes[mesh_p1.element_edges[:, k]]
        forward = t[:, i] < t[:, j]
        elements[:, 3 + 2*k:5 + 2*k] = np.where(
            forward[:, None], pair, pair[:, ::-1]
        )
    elements[:, 9] = interior_nodes

    boundary_nodes = np.unique(np.concatenate((
        mesh_p1.boundary_edges.ravel(),
        edge_nodes[mesh_p1.boundary_edge_ids].ravel(),
    )))
    arrays = dict(
        points=points, elements=elements, vertices=np.arange(nv),
        edges=edges, edge_nodes=edge_nodes, interior_nodes=interior_nodes,
        element_edges=mesh_p1.element_edges,
        boundary_edges=mesh_p1.boundary_edges, boundary_nodes=boundary_nodes,
        p1_points=p, p1_elements=t,
    )
    return P3Mesh(**{name: _readonly(value) for name, value in arrays.items()})


def plot_mesh(mesh, *, ax=None, show_nodes=False):
    """Plot affine triangles; optionally distinguish P3 node types."""
    import matplotlib.pyplot as plt

    if isinstance(mesh, P3Mesh):
        p, t = mesh.p1_points, mesh.p1_elements
    elif isinstance(mesh, P1Mesh):
        p, t = mesh.points, mesh.elements
    else:
        raise TypeError("mesh must be a P1Mesh or P3Mesh.")
    if ax is None:
        _, ax = plt.subplots()

    ax.triplot(p[:, 0], p[:, 1], t, color="0.65", linewidth=0.5)
    if show_nodes:
        if isinstance(mesh, P3Mesh):
            for ids, label, color in (
                (mesh.vertices, "Vertices", "tab:blue"),
                (mesh.edge_nodes.ravel(), "Edge nodes", "tab:orange"),
                (mesh.interior_nodes, "Interior nodes", "tab:green"),
            ):
                ax.scatter(mesh.points[ids, 0], mesh.points[ids, 1],
                           s=12, label=label, color=color)
            ax.legend()
        else:
            ax.scatter(p[:, 0], p[:, 1], s=12)
    ax.set(xlabel="R [m]", ylabel="Z [m]")
    ax.set_aspect("equal", adjustable="box")
    return ax
