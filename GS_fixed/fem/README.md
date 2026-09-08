# P1 FEM stage

This folder discretizes the fixed-boundary Grad-Shafranov weak form

```text
integral (1/R) grad(psi) . grad(v) dR dZ
    = integral [mu0 R p'(psi) + FF'(psi)/R] v dR dZ.
```

## Files

- `assemble_stiffness.m` assembles the sparse symmetric matrix for the
  `1/R`-weighted elliptic operator using three-point triangle quadrature.
- `assemble_rhs.m` accepts a constant, nodal data, `q(R,Z)`, or nonlinear
  `q(R,Z,psi)` callback and assembles the P1 load vector.
- `apply_dirichlet.m` applies scalar, nodal, or functional LCFS values by
  symmetric row/column elimination.
- `verify_fem.m` checks symmetry, the constant null mode, source integration,
  positive definiteness after Dirichlet elimination, and a manufactured
  `psi=R` solution on a concave C-shaped domain.

The source passed to `assemble_rhs` is the **weak-form source**

```text
q = mu0 R p'(psi) + FF'(psi)/R,
```

not the undivided strong source `mu0 R^2 p' + FF'`.

From `GS_fixed`, run:

```matlab
addpath('geometry', 'fem')
report = verify_fem(true);
```

Minimal assembly example:

```matlab
K = assemble_stiffness(mesh);
rhs = assemble_rhs(mesh, @(R,Z,psi) mu0*R.*pprime(psi) ...
                                   + FFprime(psi)./R, psiCurrent);
[Kbc, rhsBc, bc] = apply_dirichlet(mesh, K, rhs, 0.0);
psiNew = Kbc \ rhsBc;
```
