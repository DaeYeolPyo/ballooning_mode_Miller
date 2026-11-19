import numpy as np

from gs_fixed.fem import P3_REFERENCE_NODES
from gs_fixed.mesh import (
    P1Mesh,
    build_p1_mesh,
    generate_p3_nodes,
)


def check_p3_geometry(p1, p3):
    expected_count = (
        p1.n_vertices
        + 2 * len(p1.edges)
        + p1.n_elements
    )

    assert p3.n_dof == expected_count

    np.testing.assert_array_equal(
        p3.elements[:, :3],
        p1.elements,
    )

    xi = P3_REFERENCE_NODES[:, 0]
    eta = P3_REFERENCE_NODES[:, 1]

    barycentric = np.column_stack((
        1 - xi - eta,
        xi,
        eta,
    ))

    expected = np.einsum(
        "ij,ejk->eik",
        barycentric,
        p1.points[p1.elements],
    )

    np.testing.assert_allclose(
        p3.points[p3.elements],
        expected,
        rtol=0,
        atol=1e-13,
    )

    assert not np.intersect1d(
        p3.boundary_nodes,
        p3.interior_nodes,
    ).size


def test_two_triangles():
    points = np.array([
        [2.0, 0.0],
        [3.0, 0.0],
        [3.0, 1.0],
        [2.0, 1.0],
    ])

    # First triangle is deliberately clockwise.
    elements = np.array([
        [0, 2, 1],
        [0, 2, 3],
    ])

    p1 = P1Mesh(points, elements)
    p3 = generate_p3_nodes(p1)

    check_p3_geometry(p1, p3)

    assert p3.n_dof == 16
    assert len(p3.boundary_nodes) == 12
    assert len(p3.free_nodes) == 4

    # Shared edge: local 3->1 in element 0,
    #              local 1->2 in element 1.
    np.testing.assert_array_equal(
        p3.elements[0, 7:9],
        p3.elements[1, 3:5][::-1],
    )


def test_concave_boundary():
    # Rectangle of area 9 with a notch of area 2.
    boundary = np.array([
        [2.0, 0.0],
        [5.0, 0.0],
        [5.0, 1.0],
        [3.0, 1.0],
        [3.0, 2.0],
        [5.0, 2.0],
        [5.0, 3.0],
        [2.0, 3.0],
    ])

    p1 = build_p1_mesh(
        boundary,
        nr=13,
        nz=15,
    )
    p3 = generate_p3_nodes(p1)

    check_p3_geometry(p1, p3)

    v = p1.points[p1.elements]
    a = v[:, 1] - v[:, 0]
    b = v[:, 2] - v[:, 0]

    signed_double_area = (
        a[:, 0] * b[:, 1]
        - a[:, 1] * b[:, 0]
    )

    assert np.all(signed_double_area > 0)

    np.testing.assert_allclose(
        0.5 * signed_double_area.sum(),
        7.0,
        rtol=0,
        atol=1e-12,
    )


if __name__ == "__main__":
    test_two_triangles()
    test_concave_boundary()
    print("Mesh checks passed.")
