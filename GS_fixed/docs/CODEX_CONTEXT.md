# FEM Grad–Shafranov Solver — Codex Context

Last consolidated: 2026-09-04

## 1. Purpose and evidence level

This document records the current numerical design and implementation plan for a MATLAB fixed-boundary Grad–Shafranov (GS) solver. The immediate target is a conforming cubic (`P3`) Lagrange finite-element method on an unstructured triangular mesh in the poloidal `(R,Z)` plane.

The equilibrium solver is intended to become the front end of a longer pipeline:

```text
fixed-boundary GS equilibrium
    -> magnetic field and flux-surface geometry
    -> flux/straight-field-line coordinates and equilibrium derivatives
    -> linear ideal-MHD / ballooning stability calculations
```

Strongly shaped, indented, and C-shaped tokamak boundaries are important target cases. Accurate derivatives of `psi` will matter downstream, but the present GS discretization remains a standard `C0`, `H1`-conforming P3 method. A `C1` element is not required by the GS weak form and is not the current implementation choice.

### Important status note

At the time this file was created, the current ChatGPT project mirror contained no MATLAB source files under `sources/`. The components below were designed, reviewed, or drafted in related conversations, but their presence in the actual MATLAB repository and their runtime test results could not be verified here.

Use these labels throughout this document:

- **DECIDED**: numerical convention agreed in the conversations.
- **DRAFTED**: MATLAB implementation was written in a conversation, but is not present in this project mirror.
- **REVIEWED**: formulas/code were inspected in a conversation; this is not a substitute for executing tests.
- **VERIFIED**: reserved for a test actually run against the repository. No P3 component currently has this status in this mirror.

Before continuing implementation, locate or copy the real `.m` files into the working repository and reconcile them with the canonical interfaces below.

## 2. Governing equation and fixed-boundary weak form

Use the axisymmetric Grad–Shafranov sign convention

```text
-Delta* psi = mu0 R^2 p'(psi) + F(psi) F'(psi),

Delta* psi = R d/dR[(1/R) dpsi/dR] + d^2 psi/dZ^2.
```

Here `FFprime(psi)` means `F(psi) * dF/dpsi`, not `d(F^2)/dpsi` (the latter is twice as large).

In the two-dimensional `(R,Z)` plane,

```text
-div_RZ[(1/R) grad_RZ psi]
    = mu0 R p'(psi) + FFprime(psi)/R
    = S(R,psi).
```

For prescribed boundary flux `psi = psi_b` on `dOmega`, use trial functions satisfying the prescribed trace and test functions `v = 0` on `dOmega`. Integration by parts gives

```text
Integral_Omega (1/R) grad(v) . grad(psi) dR dZ
    = Integral_Omega v S(R,psi) dR dZ.
```

The discrete nonlinear problem is therefore

```text
K Psi = f(Psi),

K_ij = Integral_Omega (1/R) grad(N_i) . grad(N_j) dA,

f_i(Psi) = Integral_Omega N_i
             [mu0 R p'(psi_h) + FFprime(psi_h)/R] dA.
```

Consequences:

- On a fixed mesh, `K` is independent of `psi`; the nonlinearity is in the source vector.
- Before essential boundary conditions, the constant vector is in the stiffness nullspace and `K` is singular. This is expected.
- With a nonempty Dirichlet boundary and a valid domain with `R > 0`, the reduced matrix `K_II` should be symmetric positive definite.
- Keep the above sign convention consistent in profiles, manufactured solutions, residuals, and plots.

## 3. Current spatial discretization

Status: **DECIDED**

- Mesh source: a two-dimensional `delaunayTriangulation`/triangulation with three P1 vertices per element.
- Geometry: straight-sided affine triangle, even though the solution is P3.
- Solution/test space: continuous cubic Lagrange triangles with 10 local DOFs.
- Mapping: each physical element is determined only by its three vertices.
- Assembly: sparse global matrix using element triplets.
- Integration: configurable Dunavant triangle quadrature.

This is a subparametric setup: affine P1 geometry plus P3 field interpolation. The additional edge and interior P3 nodes are solution nodes; they do not make the physical boundary cubic. Curved P3 isoparametric geometry is a later extension.

## 4. Reference triangle and barycentric coordinates

Use

```text
T_hat = {(xi,eta): xi >= 0, eta >= 0, xi + eta <= 1},

lambda1 = 1 - xi - eta,
lambda2 = xi,
lambda3 = eta.
```

The physical affine map is

```text
[R; Z] = [R1; Z1] + J [xi; eta],

J = [R2-R1, R3-R1;
     Z2-Z1, Z3-Z1].
```

Enforce counter-clockwise vertex orientation and require `detJ > 0`. Do not hide an inverted element by replacing `detJ` with `abs(detJ)`. Use a scale-aware degeneracy tolerance rather than an exact `detJ == 0` comparison.

Reference-to-physical gradients satisfy

```text
grad_RZ(N) = J^(-T) grad_ref(N).
```

In MATLAB, prefer

```matlab
dNphys = J.' \ dNref;
```

to forming a general inverse. For the vectorized affine-element kernel, an explicitly derived 2-by-2 inverse transpose is acceptable if it is tested carefully.

## 5. Canonical P3 local node ordering

Status: **DECIDED; do not change silently**

| Local node | `(lambda1,lambda2,lambda3)` | `(xi,eta)` | Meaning |
|---:|---:|---:|---|
| 1 | `(1,0,0)` | `(0,0)` | vertex 1 |
| 2 | `(0,1,0)` | `(1,0)` | vertex 2 |
| 3 | `(0,0,1)` | `(0,1)` | vertex 3 |
| 4 | `(2/3,1/3,0)` | `(1/3,0)` | first node on directed edge `1 -> 2` |
| 5 | `(1/3,2/3,0)` | `(2/3,0)` | second node on directed edge `1 -> 2` |
| 6 | `(0,2/3,1/3)` | `(2/3,1/3)` | first node on directed edge `2 -> 3` |
| 7 | `(0,1/3,2/3)` | `(1/3,2/3)` | second node on directed edge `2 -> 3` |
| 8 | `(1/3,0,2/3)` | `(0,2/3)` | first node on directed edge `3 -> 1` |
| 9 | `(2/3,0,1/3)` | `(0,1/3)` | second node on directed edge `3 -> 1` |
| 10 | `(1/3,1/3,1/3)` | `(1/3,1/3)` | element interior node |

Compactly,

```text
[1 2 3 4 5 6 7 8 9 10]
 = [v1 v2 v3 e12a e12b e23a e23b e31a e31b c].
```

The `3 -> 1` direction for nodes 8 and 9 is a common source of errors. Two neighboring elements traverse a shared edge in opposite local directions; they must use the same two global edge-node IDs in reversed local order.

## 6. P3 shape-function convention

Status: **DRAFTED and mathematically REVIEWED; runtime verification pending**

Define

```text
B1(x)     = 1/2 x (3x-1)(3x-2),
B2(x,y)   = 9/2 x y (3x-1),
B3(x,y,z) = 27 x y z.
```

The canonical basis order is

```text
N1  = B1(lambda1)
N2  = B1(lambda2)
N3  = B1(lambda3)
N4  = B2(lambda1,lambda2)
N5  = B2(lambda2,lambda1)
N6  = B2(lambda2,lambda3)
N7  = B2(lambda3,lambda2)
N8  = B2(lambda3,lambda1)
N9  = B2(lambda1,lambda3)
N10 = B3(lambda1,lambda2,lambda3)
```

`construct_P3_shape_functions(xi,eta)` should return

```text
N      : 1 x 10
dNref  : 2 x 10
dNref(1,:) = dN/dxi
dNref(2,:) = dN/deta
```

Known reviewed bug: a draft ended with

```matlab
dNref = [dN_dxi.'; dN_deta.'];
```

which produces `20 x 1`. The required line is

```matlab
dNref = [dN_dxi; dN_deta];
```

The analytical derivative formulas were reviewed as consistent with the basis ordering, but they still require finite-difference and polynomial-reproduction tests in the actual repository.

## 7. `upgradeP1toP3` mesh upgrade

Status: **DRAFTED; verification code DRAFTED; not executed in this mirror**

Canonical interface:

```matlab
meshP3 = upgradeP1toP3(DT)
```

Expected output fields:

| Field | Meaning |
|---|---|
| `points` | all P3 solution-node coordinates, `Ndof x 2` |
| `elements` | P3 connectivity, `Nt x 10`, in the canonical local order |
| `vertices` | original P1 vertex IDs |
| `edges` | unique canonical P1 edges, stored low-ID to high-ID |
| `edgeNodes` | two global P3 node IDs per unique edge, following `edges` orientation |
| `interiorNodes` | one global node ID per triangle |
| `elementEdges` | global edge IDs for local edges `[1-2,2-3,3-1]` |
| `boundaryEdges` | P1 free-boundary edges |
| `boundaryNodes` | boundary vertices plus both P3 nodes on every boundary edge |
| `P1points` | original geometric vertices |
| `P1elements` | counter-clockwise P1 connectivity used by the affine map |

Global numbering is intended to be

```text
1 ... Nv                  original vertices
Nv+1 ... Nv+Ne            first node on each canonical edge
Nv+Ne+1 ... Nv+2Ne        second node on each canonical edge
Nv+2Ne+1 ... Nv+2Ne+Nt    one interior node per element
```

Therefore the invariant node count is

```text
Ndof = Nv + 2 Ne + Nt.
```

For a canonical edge `(i,j)` with `i < j`, its nodes are located at

```text
(2/3) x_i + (1/3) x_j,
(1/3) x_i + (2/3) x_j.
```

Element connectivity must reverse these two global IDs when the local directed edge is opposite to the canonical edge direction. Each element interior node is its vertex centroid.

### Mesh-upgrade verification

The proposed `validateUpgradeP1toP3` test should check all of the following:

1. `size(elements) == [Nt,10]` and all IDs are valid integers.
2. `Ndof == Nv + 2*Ne + Nt`.
3. The first three local IDs match the counter-clockwise P1 element vertices.
4. Every element has positive, nondegenerate `detJ`.
5. Local nodes 4–9 are at the exact directed one-third/two-third edge positions.
6. Local node 10 is the element centroid.
7. Both elements adjacent to an interior edge share exactly the same two global edge-node IDs, in opposite local order when appropriate.
8. Each element has one unique interior node; no interior node is shared.
9. Boundary-node classification equals boundary vertices plus two nodes per unique boundary edge.
10. There are no accidental duplicate coordinates/IDs beyond intended shared-edge nodes.

Visual inspection should use a small connected one-ring patch rather than plotting every label on a large mesh. Show vertex, edge, and interior nodes with different markers, and display both local and global IDs on the seed element and its neighbors. This specifically exposes edge-orientation mistakes.

## 8. Dunavant quadrature

Status: **rules for degrees 1–20 DRAFTED; implementation review found a permutation bug; verifier DRAFTED; no recorded execution result**

Canonical interface:

```matlab
[xi_q, eta_q, w_q] = Gauss_Dunavant_quadrature(p)
```

The discussed implementation expands symmetric barycentric orbits:

- centroid: one point;
- type 1: permutations of `(a,b,b)`, three points;
- type 2: permutations of `(a,b,c)`, six points.

With `(lambda1,lambda2,lambda3)` mapped to `(xi,eta) = (lambda2,lambda3)`, the correct type-1 helper is

```matlab
xi  = [b, a, b];
eta = [b, b, a];
```

A draft used `xi=[a,b,a]`, `eta=[b,a,a]`; this is wrong for the documented `(a,b,b)` orbit and should fail the degree-2 moment test. The discussed type-2 permutation helper was reviewed as consistent.

### Weight convention — resolve before assembly

Two conventions appeared in the conversations:

- Dunavant table weights normalized so `sum(w) = 1`;
- reference-triangle integration weights scaled so `sum(w) = 1/2`.

The affine element formulas in this document assume direct reference-triangle weights:

```text
sum(w) = area(T_hat) = 1/2,
dA = detJ dxi deta.
```

If `Gauss_Dunavant_quadrature` returns normalized weights with `sum(w)=1`, multiply them by `1/2` exactly once when constructing the canonical quadrature object. Never compensate in some element routines but not others.

### Exactness verification

`verify_Gauss_Dunavant_quadrature` was drafted to test `p=1,...,20` using

```text
Integral_T xi^i eta^j dxi deta = i! j! / (i+j+2)!,
for every i,j >= 0 with i+j <= p.
```

It should also check:

- expected point counts: `[1,3,4,6,7,12,13,16,19,25,27,33,37,42,48,52,61,70,73,79]`;
- weight normalization;
- equal lengths and finite values of `xi`, `eta`, and `w`;
- barycentric-coordinate consistency;
- number of negative weights and exterior points, reported but not automatically treated as errors for high-degree Dunavant rules;
- exactness at degree `p`, and optionally failure/non-guarantee at degree `p+1`.

Do not infer that an exterior point or negative weight alone makes a high-order rule invalid.

### Quadrature degree policy

An early stiffness draft used the 7-point degree-5 rule. For P3,

```text
grad(N_i) . grad(N_j) is degree 4,
N_i N_j is degree 6,
```

but the GS stiffness integrand also contains `1/R`, and nonlinear/curved-geometry terms are not polynomials. Use a configurable Dunavant degree; degree 7 is a sensible starting point. Establish quadrature convergence by comparing at least degrees 7, 8, and 10. Do not confuse “7 points” with “degree 7 exactness”: the 7-point rule is degree 5.

## 9. Precomputed quadrature object

Status: **DRAFTED design**

Create the basis once per selected quadrature rule:

```text
quad.xi      : nq x 1
quad.eta     : nq x 1
quad.w       : nq x 1, with sum(w)=1/2
quad.N       : nq x 10
quad.dNdxi   : nq x 10
quad.dNdeta  : nq x 10
```

Avoid maintaining incompatible alternatives such as both a `2 x 10 x nq` `quad.dNref` and separate derivative matrices unless there is a clear adapter. The preferred vectorized affine kernel uses `dNdxi` and `dNdeta` as `nq x 10` arrays.

## 10. Element stiffness and global sparse assembly

Status: **local and global routines DRAFTED; vectorized kernel DRAFTED; not executed in this mirror**

For one affine triangle,

```text
Ke = Sum_q w_q detJ/R_q [grad(N)_q]^T [grad(N)_q],

R_q = lambda1_q R1 + lambda2_q R2 + lambda3_q R3.
```

With physical derivative matrices `dNdR,dNdZ` of size `nq x 10` and

```matlab
alpha = quad.w .* detJ ./ Rq;
```

the preferred element kernel is

```matlab
Ke = dNdR.' * (alpha .* dNdR) ...
   + dNdZ.' * (alpha .* dNdZ);
```

Optionally remove roundoff asymmetry with `Ke = 0.5*(Ke+Ke.')`, but first test the uncleaned matrix during debugging so symmetry errors are not hidden.

Global assembly should keep the element loop and preallocate `100*Nt` triplets:

```text
I, J, V -> sparse(I,J,V,Ndof,Ndof).
```

Each element uses `meshP3.P1points(meshP3.P1elements(e,:),:)` for affine geometry and `meshP3.elements(e,:)` for its 10 solution DOFs.

## 11. Element source vector and global RHS

Status: **DRAFTED; not executed in this mirror**

Canonical local interface:

```matlab
fe = construct_P3_source_vector( ...
    Xe, psi_e, quad, pprime_func, FFprime_func)
```

At quadrature points,

```text
psi_q = quad.N * psi_e,
R_q   = lambda1_q R1 + lambda2_q R2 + lambda3_q R3,
S_q   = mu0 R_q pprime_func(psi_q)
      + FFprime_func(psi_q)/R_q,

fe = quad.N.' * ((quad.w .* detJ) .* S_q).
```

Require profile function handles to accept and return vector-shaped inputs. Convert returned profiles to columns and validate their lengths. Check `R_q > 0` and finite source values.

Canonical global interface:

```matlab
f = assemble_P3_source_vector( ...
    meshP3, psi, quad, pprime_func, FFprime_func)
```

Direct accumulation

```matlab
f(ids) = f(ids) + fe;
```

is acceptable initially. An `accumarray` version with `10*Nt` entries was also drafted for larger meshes. Since the source depends on `psi`, global RHS assembly is repeated at every nonlinear iteration.

For a constant source `S0`, the essential element test is

```text
sum(fe) = S0 * area(element) = S0 * detJ/2.
```

This test simultaneously catches partition-of-unity, Jacobian, and quadrature-weight scaling errors.

## 12. Dirichlet boundary condition and solves

Status: **reduction routine DRAFTED; complete P3 solve not yet verified**

For boundary IDs `B` and interior IDs `I`, prescribe `Psi_B` and solve

```text
K_II Psi_I = f_I - K_IB Psi_B.
```

For homogeneous fixed boundary, `Psi_B=0` and the correction vanishes. Restore a full global vector after solving. Validate that the ordering of a vector-valued `Psi_B` matches `boundaryNodes` exactly.

### First linear solve

Start with constant profiles,

```matlab
pprime_func  = @(psi) Cp * ones(size(psi));
FFprime_func = @(psi) CF * ones(size(psi));
```

so that the RHS is independent of `psi`. Complete a manufactured-solution or benchmark solve before attempting nonlinear profiles.

### Nonlinear solve

Use Picard iteration first:

```text
K Psi_solve^(k+1) = f(Psi^k),

Psi^(k+1) = (1-omega) Psi^k + omega Psi_solve^(k+1),
0 < omega <= 1.
```

Factorize the fixed reduced stiffness matrix once and reuse the factorization. Monitor both a solution-update norm and the true nonlinear residual

```text
r(Psi) = K Psi - f(Psi)
```

on the free DOFs. Add iteration limits, stagnation/divergence detection, and under-relaxation.

Newton is a later option. Its tangent contains

```text
K_ij - Integral N_i [dS/dpsi] N_j dA,

dS/dpsi = mu0 R p''(psi) + d(FFprime)/dpsi / R.
```

Do not add Newton until the linear and Picard paths have independent tests.

## 13. Vectorization policy

Status: **DECIDED**

Prioritize numerical transparency and memory safety:

1. Precompute `N`, `dNdxi`, and `dNdeta` at quadrature points once.
2. Vectorize the quadrature direction inside each affine element using matrix products.
3. Keep the global element loop initially.
4. Use preallocated sparse triplets for matrices.
5. Use direct RHS accumulation first; switch to `accumarray` only after profiling.
6. Do not create dense global matrices or a dense diagonal quadrature-weight matrix.
7. Avoid whole-mesh vectorization if it obscures local-to-global indexing or creates large temporary 3-D arrays.
8. Optimize only after the verified scalar/loop implementation and vectorized implementation agree to tolerance.

## 14. Canonical MATLAB files and roles

The following files were named or drafted in the conversations. Confirm which ones exist in the real repository before editing.

| File | Role | Current evidence |
|---|---|---|
| `upgradeP1toP3.m` | convert P1 triangulation to global 10-node P3 connectivity and boundary-node sets | DRAFTED |
| `construct_P3_shape_functions.m` | evaluate canonical P3 basis and reference derivatives | DRAFTED, formula REVIEWED; known output-shape fix |
| `Gauss_Dunavant_quadrature.m` | return Dunavant rules for requested degree 1–20 | DRAFTED; known type-1 orbit fix |
| `construct_P3_quadrature.m` | normalize weights and precompute P3 basis/derivative arrays | DRAFTED design |
| `construct_P3_stiffness_matrix.m` | assemble one `10 x 10` affine-element GS stiffness matrix | DRAFTED |
| `assemble_P3_stiffness_matrix.m` | assemble global sparse stiffness using triplets | DRAFTED |
| `construct_P3_source_vector.m` | assemble one nonlinear GS element load vector | DRAFTED |
| `assemble_P3_source_vector.m` | assemble the global nonlinear RHS | DRAFTED |
| `apply_P3_dirichlet.m` | form the reduced Dirichlet system and corrected RHS | DRAFTED |
| `verify_Gauss_Dunavant_quadrature.m` | verify rules by all monomial moments through degree `p` | DRAFTED |
| `validateUpgradeP1toP3.m` | verify counts, positions, topology, sharing, and boundary classification | DRAFTED |
| `visualizeP3Mesh.m` | plot P3 node types and global numbering on a small mesh | DRAFTED helper |
| `visualizeP3Element.m` | inspect local nodes 1–10 on one element | DRAFTED helper |
| `visualizeP3Patch.m` | inspect a connected one-/two-ring patch and shared-edge IDs | DRAFTED helper |

Names such as `shapeP3.m` and `triangle_gauss_points.m` appeared as earlier alternatives. Prefer the canonical names above and avoid duplicate implementations with different conventions. `GS_source.m` may be introduced to isolate the physical profile model from FEM integration once the core source-vector test passes.

## 15. Numerical conventions

- Coordinates are `(R,Z)` in that order; require `R > 0` at all relevant nodes and quadrature points.
- Element vertices are counter-clockwise; `detJ > 0` is an invariant.
- Reference triangle area is `1/2`.
- Canonical quadrature object weights sum to `1/2`.
- `N` is `nq x 10` in precomputed storage; a single-point call may return `1 x 10`.
- `dNref` from the shape function is `2 x 10`; precomputed `dNdxi,dNdeta` are `nq x 10`.
- Global nodal vectors are columns.
- Geometry remains affine until a separate curved-geometry design is implemented.
- Use the local node ordering in Section 5 in every mesh, basis, assembly, interpolation, and plotting routine.
- `FFprime` means `F*dF/dpsi`.
- Keep SI units and `mu0 = 4*pi*1e-7` explicit unless a documented nondimensionalization is introduced.
- Tolerances should be relative to an appropriate matrix/vector/geometry scale, not only absolute constants.

## 16. Verification checklist

Do not call the P3 solver complete until these checks are automated and passing.

### Shape basis

- [ ] Kronecker property `N_i(x_j)=delta_ij` at all 10 reference nodes.
- [ ] Partition of unity: `sum(N,2)=1` at vertices, edges, interior random points, and quadrature points.
- [ ] Derivative sum: `sum(dNref,2)=[0;0]`.
- [ ] Central finite-difference derivatives agree with analytical `dNref`.
- [ ] Reproduce every polynomial of total degree `<=3` from nodal interpolation.

### Mesh upgrade

- [ ] Node-count and connectivity-size identities.
- [ ] Positive orientation and no degenerate elements.
- [ ] All edge and centroid positions match the canonical barycentric coordinates.
- [ ] Interior shared edges reuse the same global node pair with correct reversal.
- [ ] Boundary node set is complete and contains no element-interior nodes.
- [ ] One-ring visual checks at interior, boundary, and highly distorted regions.

### Quadrature

- [ ] Expected point count for degrees 1–20.
- [ ] Canonical weight sum equals `1/2` after normalization.
- [ ] All monomial moments through requested degree pass.
- [ ] Exterior points/negative weights are reported and understood.
- [ ] Solver results are insensitive to increasing degree from 7 to 8/10 within target tolerance.

### Element and global stiffness

- [ ] `Ke` and `K` are symmetric before optional symmetry cleanup.
- [ ] `Ke*ones(10,1)` and `K*ones(Ndof,1)` are near zero before Dirichlet reduction.
- [ ] `Ke` is positive semidefinite with one constant null mode.
- [ ] `K_II` is positive definite after Dirichlet reduction (e.g. Cholesky succeeds).
- [ ] Triplet assembly matches a simple reference assembly on a tiny mesh.
- [ ] Affine-coordinate scaling and rotated/skewed triangle tests pass.

### Source and boundary condition

- [ ] Constant-source element identity `sum(fe)=S0*area`.
- [ ] Direct and `accumarray` global RHS assemblies agree.
- [ ] Constant-profile source is independent of the supplied `psi` vector.
- [ ] Homogeneous and nonzero Dirichlet reductions match a direct row/column reference treatment.
- [ ] Nonlinear residual is evaluated only after restoring prescribed boundary values.

### End-to-end solution

- [ ] Manufactured solution on a domain bounded away from `R=0`.
- [ ] `h`-refinement convergence study for `psi` and `grad(psi)`.
- [ ] `p`/quadrature convergence checks.
- [ ] Flux and source sign conventions confirmed against an independent solution or analytic case.
- [ ] P3 visualization samples each element internally; a vertex-only `trisurf` is not accepted as a P3 field plot.

## 17. Coding conventions

- Use MATLAB functions for reusable components; keep demonstration/driver scripts thin.
- Preserve canonical function signatures and local ordering unless a coordinated migration updates every consumer and test.
- Validate dimensions, finiteness, orientation, index ranges, and profile output shapes at public interfaces.
- Prefer sparse matrices and preallocation.
- Avoid `inv` for general linear solves.
- Keep mathematical symbols recognizable in variable names (`detJ`, `dNdxi`, `dNdR`, `psi_q`, `pprime_q`).
- Document array dimensions beside nontrivial matrix operations.
- Separate mesh topology, reference basis, geometry mapping, quadrature, physics source, assembly, boundary conditions, nonlinear solve, and post-processing.
- Make tests deterministic; use a fixed random seed when sampling elements/points.
- Do not weaken a test merely to make an implementation pass. Explain any tolerance change with scale and precision evidence.
- Profile before introducing complex whole-mesh vectorization.

## 18. Known cautions and open TODOs

### Must fix or confirm first

1. Find the actual MATLAB repository/files; this project mirror currently contains only context material.
2. Confirm the `dNref` return is `2 x 10`, without transposed concatenation.
3. Correct and test the Dunavant `(a,b,b)` orbit expansion.
4. Standardize quadrature weights to `sum(w)=1/2` in `construct_P3_quadrature`.
5. Confirm every routine uses the same node ordering, especially nodes 8 and 9 on edge `3 -> 1`.
6. Replace exact zero-area comparisons with a scale-aware degeneracy check.

### Next implementation sequence

1. Put the canonical P3 source files under version control and remove/adapter-wrap duplicate historical names.
2. Run basis, mesh-upgrade, and Dunavant unit tests before assembling the PDE.
3. Implement/verify the precomputed quadrature object.
4. Verify loop and vectorized local stiffness kernels against each other.
5. Verify sparse global stiffness assembly on a two-triangle shared-edge mesh.
6. Verify local/global source assembly with `S=1` and constant `p'`, `FFprime`.
7. Apply homogeneous and nonzero Dirichlet conditions and solve a linear manufactured problem.
8. Add P3-aware sampling/contour visualization.
9. Add Picard iteration with reusable factorization, relaxation, and residual monitoring.
10. Validate a nonlinear GS benchmark; add Newton only if Picard performance requires it.
11. Run mesh and quadrature convergence on strongly shaped boundaries.
12. Add derivative recovery or projection for smooth `B`, curvature, metric, and ballooning inputs if element-boundary gradient jumps materially affect downstream results.

### Later extensions, not current scope

- Curved P3 isoparametric boundary elements. Then all 10 geometry nodes participate, `J` varies at every quadrature point, and `detJ>0` must be checked throughout each element.
- Boundary-edge projection onto an analytic or tabulated plasma boundary.
- Adaptive `h/p` refinement.
- Smooth gradient/Hessian recovery or projection to a `C1/C2` spline field for downstream equilibrium derivatives.
- Coupling the verified GS equilibrium to flux-coordinate and linear ideal-MHD/ballooning solvers.
- Structure-preserving/compatible finite-element spaces for vector MHD variables. This is a separate linear-MHD discretization question and should not be conflated with the present scalar P3 GS solver.

## 19. Recommended starting prompt for a new Codex task

```text
Read AGENTS.md and CODEX_CONTEXT.md first. Locate the actual MATLAB P3
Grad–Shafranov source files, compare them with the canonical conventions in
CODEX_CONTEXT.md, and report any mismatch before changing numerical behavior.
Start by running or completing the shape-function, mesh-upgrade, and Dunavant
quadrature verification checks. Preserve the P3 local node ordering and the
reference-triangle weight convention.
```
