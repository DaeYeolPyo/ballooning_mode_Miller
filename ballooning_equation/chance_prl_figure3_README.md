# Chance PRL (1983), Figure 3: native solver benchmark

This benchmark uses the boundary and prescribed profiles in the supplied
Chance, Jardin and Stix paper, PRL **51**, 1963–1966,
DOI `10.1103/PhysRevLett.51.1963`. It does not use the outboard indented
Miller boundary of the C-shape scan.

## Paper inputs

Let `u = psiN` increase from the magnetic axis (0) to the boundary (1).
The profile variable is normalized flux, **not** `sqrt(psiN)`.

```text
rho(t) = A [1 + B cos(t)]
gamma(t) = C sin(t)
R(t) = xbar + rho(t) cos(gamma(t))
Z(t) = E rho(t) sin(gamma(t)),   0 <= t < 2*pi

A = 1, B = 0.6, Rgeo/a = 4, b/a = 1.386
p(u) = p0 (1-u^2)^2
dp/du = -4 p0 u (1-u^2)
q(u) = 1.03 + 0.84375 u - 1.1775 u^2 + 3.50375 u^3
```

For each requested inboard indentation `d/(2a)`, solve for `C`, then adjust
`E` and `xbar` to keep the two geometric ratios fixed. Keep `A=1`; the actual
minor radius changes with `C`. The indentation depth is the difference
between `R(pi)` and the minimum boundary radius. `chance_prl_bean_boundary`
implements these definitions and analytic first/second boundary derivatives.

The native fixed-boundary GS equation uses the existing P3 FEM formulation.
Prescribe `dp/du` and `q(u)` and solve for the FFprime coefficients. In this
solver FFprime means `F dF/du`; division by the dimensional native flux span
converts it to `F dF/dpsi`. Set `Fboundary = Rgeo * BReference`. The chosen
field scale is 3 T; the dimensionless equilibrium/stability comparison is
made with the paper's volume beta:

```text
dV = 2*pi*R*dR*dZ
beta = 2*mu0*integral(p*dV) / integral((Bp^2+Bphi^2)*dV)
```

## Calculation and acceptance

1. Generate the Eq. (2) boundary and prepare the GS operators once per shape.
2. Vary physical `p0`, using converged equilibria as continuation predictors.
   Compute volume beta after every solve; central beta is not the plot ordinate.
3. Accept only equilibria meeting the GS and prescribed-q tolerances.
4. Construct each flux surface from the actual FEM field and calculate direct
   geometry, curvature, local shear, and the ballooning coefficients.
5. Solve the Dirichlet ballooning eigenproblem over radial labels and ballooning
   phases. The largest eigenvalue is positive for instability in our convention.
6. Plot the classifications at the **measured** beta. A failed surface prevents
   a stable label; any valid positive eigenvalue still proves instability.
7. Refine the equilibrium mesh, radial/phase sampling and ballooning interval
   before claiming a converged marginal curve or stability on every flux surface.

`run_chance_prl_figure3_test` targets beta to 1.5% relative accuracy, saves all
native equilibria and beta-targeting histories, and resumes matching checkpoints.
An unmatched configuration requires a different output directory. The reference
curve in the plot is an approximate visual reading of the paper, not a fit to
the computed points. It is not original numerical data.

```matlab
setup_indented_cshape_scan
addpath('verification')
geometry = test_chance_prl_geometry;
config = make_chance_prl_figure3_config( ...
    Indentation=[.11 .185 .25 .304 .325], GridSize=[48 64]);
report = run_chance_prl_figure3_test(config, ...
    fullfile(pwd,'output','chance_prl_figure3','fine_benchmark'));
```

The small 24-by-32 exploratory mesh is not an accepted resolution for stable
classification: some inner surfaces exceed the stricter 1% closure criterion.

## Diagnostic finding: constrain q beyond the tested surface

The original scan fitted q through `u=.95` and tested stability there. The
actual radial derivative near the end of that fit can disagree substantially
with the prescribed derivative, even when the q values themselves pass.

For the `d/(2a)=.304` low-beta fine-mesh test, fitting only through .95 gave
an independently measured `dq/du(.95)=9.05684`, while the paper prescribes
8.09290. Increasing only the Fourier representation did not fix the closure.
Extending the q constraints through .975 reduced the .95 surface closure
defect from 7.09% to 0.0917%, with a complete stable calculation.

The benchmark therefore fits q through .975 and samples ballooning through
.95. It retains the actual native q profile in the ballooning adapter; no
replacement by prescribed q/qprime is used to manufacture a classification.
The surface closure acceptance is 1%, tighter than the general scan's 5%.

At high pressure, Cartesian contour extraction near .975 can also fail when
the surface approaches the fixed boundary. The exploratory high-pressure
diagnostics test contour resolution and pressure-aware initial FFprime values
separately. These predictions do not bypass the final GS/q convergence checks.

The completed representative test did **not** reproduce the high-beta stable
region in Figure 3. The sampled first boundary at indentation .25 was bracketed
at beta 2.498–2.514%. At beta 14.35%, native equilibria meeting GS/q tolerances
still produced a positive eigenvalue on an outer surface, including GS 96-by-128
with 768 boundary points. Theta/Fourier refinement and a dense eigenvalue check
confirmed the numerical discrepancy. The cause needs independent validation
against paper Eq. (3) and neighboring-surface shear; it is not established here.

The benchmark defaults now use 48-by-64 GS, contour/adapter NR=800,
2048 surface points and 96 ballooning cells per turn. These are exploratory
controls, not a certified resolution for a full Figure 3 marginal curve.
Schema-1 stability checkpoints are rejected after the contour corrections.

Actual results, native equilibria, and their limitations are in
`output/chance_prl_figure3/validation/validation_report.md` and `validation.mat`.
The comparison figure also shows the coarse false-stability result and its
theta refinement. Run `build_chance_validation_report` to regenerate these
artifacts from the saved experiments under `tmp/chance_prl`.
