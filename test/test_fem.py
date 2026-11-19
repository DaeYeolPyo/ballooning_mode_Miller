import numpy as np

from gs_fixed.fem import (
    P3_REFERENCE_NODES,
    construct_p3_shape_functions,
    precompute_quadrature,
    construct_k_matrix,
)

# Kronecker property: reference nodes에서 N = I.
nodes = P3_REFERENCE_NODES
N, _ = construct_p3_shape_functions(nodes[:, 0], nodes[:, 1])
np.testing.assert_allclose(N, np.eye(10), atol=1e-14)

# 반시계 방향의 물리 삼각형 [R, Z].
vertices = np.array([
    [2.0, -0.3],
    [3.0,  0.2],
    [2.2,  0.9],
])

quad = precompute_quadrature(p=10)
Ke = construct_k_matrix(vertices, quad)

print(quad.nq)    # 25
print(Ke.shape)   # (10, 10)

np.testing.assert_allclose(quad.w.sum(), 0.5, atol=1e-13)
np.testing.assert_allclose(Ke, Ke.T, atol=1e-13)
np.testing.assert_allclose(Ke @ np.ones(10), 0.0, atol=1e-13)
