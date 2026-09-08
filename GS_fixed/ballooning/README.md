# Ballooning stage 1: local flux-surface geometry

This stage connects the cosine fixed-boundary FEM equilibrium to the future
local ballooning calculation without creating a GEQDSK file.

`extract_local_surface_bundle.m` traces several nested P1 flux contours
directly on the triangular plasma mesh. Each surface is started at the
outboard-midplane crossing, oriented upward, and sampled at uniform geometric
arclength using

```text
vartheta = 2*pi*s/L.
```

This `vartheta` is periodic but is not yet a straight-field-line coordinate.
No safety factor, magnetic shear, alpha, ballooning coefficient, or
eigenvalue is constructed in stage 1.

After running `main_cos`, verify this stage with:

```matlab
addpath('ballooning')
report = verify_local_surface_bundle(equilibrium,true);
```

The default verification extracts `psiN=0.55:0.10:0.95` with 512 points per
surface. These five surfaces provide the minimum radial neighborhood needed
by the existing five-point ballooning-coefficient machinery in a later stage.

Because the equilibrium uses P1 FEM, its raw contours are piecewise linear.
Increasing only `nPoloidal` adds points to the same polygon and does not remove
visible corners. The verified `main_cos` setting therefore uses
`cfg.mesh.targetH=0.05`; `nPoloidal=512` then provides a dense common periodic
grid without pretending to add geometric information absent from the FEM mesh.

## Stage 2: local equilibrium and PEST angle

`build_local_pest_equilibrium.m` evaluates the recovered FEM field on every
stage-1 surface, reconstructs `F(psi)` from normalized-power or tabulated
`FFprime`/`F` data, and computes

```text
q = (1/2*pi) integral F/(R |grad psi|) dl,
V' = 2*pi integral R/|grad psi| dl,
sqrt(g)_PEST = q R^2/F.
```

It calls the existing generic functions in the sibling `equilibrium` folder
without modifying that folder. The resulting surfaces are remapped to a
uniform 512-point PEST straight-field-line grid. The stage-2 verification
uses nine radial surfaces at `psiN=0.55:0.05:0.95`, so the target
`psiN=0.75` has a centered five-surface derivative stencil. Verify with:

```matlab
report = verify_local_pest_equilibrium(equilibrium,true);
```

This stage reports the equilibrium `q`, `qprime`, `sHat`, and reference
`alpha`, but does not yet construct ballooning coefficients or solve an
eigenvalue problem.

## Stage 3: PEST metric and ballooning coefficients

`build_local_ballooning_coefficients.m` converts the SI FEM equilibrium to
the dimensionless convention used by the existing sibling routines:

```text
Rbar   = R/aN,              aN = magnetic-axis R,
psibar = psi/(BN*aN^2),     BN = |F_axis|/aN,
pbar   = mu0*p/BN^2.
```

It then calls `equilibrium/compute_metrics.m` and
`ballooning_equation/calculate_ballooning_coeffs.m` without modifying either
sibling folder. The resulting arrays satisfy

```text
d/dtheta (g dX/dtheta) + c X = lambda f X.
```

The P1 contour and recovered FEM field contain small mesh-scale ripples that
are strongly amplified by the secular `qprime*(theta-theta0)` terms. Before
forming the metric, the stage therefore retains Fourier modes `|m|<=24` in
the PEST `R`, `Z`, `lambda`, and `B^2` arrays. The cutoff is exposed as
`poloidalModeCutoff`; verification also limits the resulting surface
displacement and the RMS change in `B^2` so this cannot silently reshape the
equilibrium.

Verify the metric identities, SI-to-dimensionless conversion, coefficient
positivity, and the target `psiN=0.75` surface with:

```matlab
report = verify_local_ballooning_coefficients(equilibrium,true);
```

This stage intentionally stops before matrix assembly or an eigenvalue
solve. `nPeriods` in `build_local_ballooning_coefficients` is already exposed
so that a later stage can construct a sufficiently long field-line domain.

## Stage 4: local ballooning eigenmode

`solve_local_ballooning_mode.m` approximates the infinite field line by
`[-nPeriods*pi,nPeriods*pi]`, imposes `X=0` at both ends, and reuses the P2
finite-element assembler in `sturm_liouville_solver/construct_matrix.m`.
Because the existing generic eigenvalue function converts sparse matrices to
dense form, the GS wrapper instead uses a sparse symmetric generalized
`eigs` solve for only the leading modes. With the adopted equation convention,
`lambda>0` is unstable and `lambda=0` is marginal.

Run the complete stage verification with:

```matlab
report = verify_local_ballooning_eigensolver(equilibrium,true);
```

The verification compares the sparse result with the existing dense solver
on a small identical matrix, then checks convergence with both field-line
length and angular resolution. This stage evaluates only the equilibrium
point at `psiN=0.75`; it does not yet vary independent `s` and `alpha`.

There are two distinct field-line-length checks. A positive isolated
eigenvalue must converge and its mode must decay before the Dirichlet ends.
At a stable point there is instead no positive discrete mode: the largest
finite-domain value approaches the continuum edge as
`lambda1 ~ -constant/nPeriods^2`. The verifier detects and reports these two
cases separately rather than incorrectly demanding a localized eigenfunction
from a stable continuum state.

## Stage 5: frozen-geometry s-alpha diagram

`prepare_local_salpha_model.m` freezes the verified PEST geometry and maps
independent scan variables through the repository contour definitions,

```text
qprime = sHat*q*Vprime/(2*V),
pprime = -alpha*(2*pi^2/Vprime)*sqrt(2*pi^2*Raxis/V).
```

Thus the equilibrium's own `(sHat,alpha)` point exactly recovers its original
`qprime`, `pprime`, and `g,c,f`. This is explicitly a frozen-geometry local
scan: `B2_psi` and the PEST metric remain those of the GS solution while
`qprime` and `pprime` vary independently. It is not a family of globally
self-consistent Grad-Shafranov equilibria far from the reference point.

The coefficient dependence is polynomial: `g,f` are quadratic in `sHat`,
while `c` contains `alpha`, `alpha^2`, and `sHat*alpha`. The scanner assembles
these sparse P2 matrix bases once per ballooning phase and then maximizes the
leading eigenvalue over `theta0`. Run the verified 17-by-17 diagram with:

```matlab
report = verify_local_salpha_scan(equilibrium,true);
```

The verification uses a regular 17-by-17 base grid and inserts the exact GS
reference `(sHat,alpha)` as an additional row and column.

The black `lambda=0` contour is the approximate local marginal boundary;
positive values are unstable. The red star marks the actual GS equilibrium
at `psiN=0.75`.

## Miller D-shape cross-check

`main_miller.m` solves a second fixed-boundary GS equilibrium using

```text
R = R0 + a*cos(theta + asin(delta)*sin(theta)),
Z = kappa*a*sin(theta),
```

with `R0/a=3.17`, `kappa=1.66`, and `delta=0.416`. It deliberately uses the
same P1 mesh target, pressure derivative, `FFprime`, and Picard settings as
`main_cos`. The complete independent cross-check is:

```matlab
run_miller_salpha_validation
```

This runs the Miller-boundary equilibrium, verifies one O-point/no X-point,
constructs the local PEST map at `psiN=0.75`, and applies exactly the same
frozen-geometry s-alpha scanner. No function in the sibling `Miller` folder
is required or modified.
