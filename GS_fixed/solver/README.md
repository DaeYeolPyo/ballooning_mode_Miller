# Solver stage

This folder solves the nonlinear fixed-boundary Grad-Shafranov problem and
tracks solutions as C-shaping is increased.

## Files

- `initial_guess.m`: solves a unit-source elliptic problem on the actual LCFS
  and scales it to a requested axis value. It fixes the LCFS exactly without
  prescribing internal flux-surface shapes.
- `solve_picard.m`: nonlinear source evaluation, linear FEM solve,
  under-relaxation, dynamic `psiAxis` update, and update/residual convergence
  tests.
- `transfer_solution.m`: P1 barycentric transfer between continuation meshes,
  with nearest-node extrapolation only for newly added plasma regions.
- `solve_continuation.m`: interpolates `A,B,C,G,H` for harmonic boundaries,
  `A,B,D,beta,G,H` for localized boundaries, or
  `A,Rback,Rtip,Rend,Zend,pInner,pOuter` for pointed boomerang boundaries;
  third-harmonic cosine boundaries interpolate `R0,A,B,C,kappa`. It remeshes
  each LCFS and follows the converged equilibrium branch.
- `verify_solver.m`: verifies a direct linear solution, a nonlinear normalized
  profile, and three-step C-shaping continuation.

Run from `GS_fixed`:

```matlab
addpath('geometry', 'fem', 'physics', 'solver')
report = verify_solver(true);
```

Basic nonlinear solve:

```matlab
K = assemble_stiffness(mesh);
psi0 = initial_guess(mesh, K, ...
    struct('boundaryValue',0,'axisValue',-1));

opts = struct('omega',0.5, 'tolerance',1e-8, ...
              'residualTolerance',1e-8, 'axisMode','min');
[psi, result] = solve_picard(mesh, profiles, psi0, opts);
```

For a negative source with `psi_boundary=0`, use a negative initial axis and
`axisMode='min'`. For the opposite flux orientation, use a positive initial
axis and `axisMode='max'`. `axisMode='auto'` selects the interior extremum
farthest from the boundary value.
