# Dense original C-shape beta map

The dense scan uses the original **outboard indented C-shape**, with fixed
measured dimensions `R_g/a_g = 4.8300`, `a_g = 0.5744 m`,
`kappa_g = 4.1781`, `delta = -0.4`, and `B_ref = 3 T`.
The dimensions are derived from the existing reference shape, rather than
changing the aspect ratio as indentation changes. Exact values are saved
in `scan.config`.

## Completed calculation (2026-10-02)

All 374 requested equilibria and their 11-by-8 stability samples completed:
147 stable, 227 unstable, zero marginal or unknown points. Actual volume
beta ranges from 0.065865% to 13.258216%. All 17 first-instability brackets
were refined with 85 additional native equilibrium/stability solves.
The first critical beta bracket decreases from 3.372583--3.385122% at
`I=0.6` to 2.467453--2.478669% at `I=1.0`. No sampled return to stability
was observed in the scanned pressure range.

The saved final-data audit checks native pressure/shape identity, GS and q
acceptance, all 88 finite phase eigenvalues per point, surface quality,
and the sign and pressure width of every transition bracket. Its result
is `output/indented_cshape_dense/quality_verification.mat`. The maximum
accepted GS residual is `4.76892e-7`; maximum constraint relative q error
is `0.00199223`, both within the unchanged limits.

The construction curve is

```text
R(theta) = R_c + a_0 [cos(theta + asin(delta) sin(theta))
                     - I ((1 + cos(theta))/2)^3],
Z(theta) = kappa_0 a_0 sin(theta).
```

The existing normalization selects `R_c`, `a_0`, and `kappa_0` to preserve
the measured box dimensions. The default horizontal axis is the measured
outboard indentation `i_out = d_out/(2 a_g)`; an additional plot uses the
construction parameter `I`.

## Sampling and equations

- `I = 0.6:0.025:1.0`: 17 boundaries.
- Pressure controls `beta_0 = [0.002, 0.01:0.01:0.12, 0.14:0.02:0.30]`:
  22 pressures per boundary, 374 requested equilibria.
- `p_0 = beta_0 B_ref^2/(2 mu_0)` and `p(u) = p_0 (1-u^2)^2`,
  where `u = psi_N`.
- Fixed target `q(u) = 1.03 + 0.84375 u - 1.1775 u^2 + 3.50375 u^3`.
- P3 GS mesh: `NR=48`, `NZ=64`, 128 boundary vertices, Dunavant degree 7.
  GS residual tolerance `5e-7`; maximum relative q error `2e-3`,
  RMS logarithmic q error `1e-3`. Failed equilibria are retained as unknown.
- The dense scan uses 12 FFprime basis functions. An 8-basis preliminary
  scan stalled near the 0.2% q tolerance for some shapes; at `I=0.65`,
  `beta_0=0.05`, 12 bases reduced maximum q error to 0.02498% in two outer
  iterations without changing the tolerances. Increasing only the q
  contour grid from 200 to 400 did not remove the 8-basis residual floor.
- Each equilibrium samples 11 radial surfaces
  `[0.05:0.1:0.85, 0.90, 0.95]` and 8 ballooning phases. The angular domain
  is `[-5*pi,5*pi]`, with 96 cells/turn, 2048 surface samples, 48 Fourier
  modes, and a 513-point radial contour-discovery grid.
- Each plotted beta is computed directly using degree-10 volume quadrature:
  `beta_av = 2 mu_0 integral(p dV)/integral((Bp^2+Bphi^2) dV)`,
  `dV = 2 pi R dR dZ`. The pressure control `beta_0` is not plotted as beta.

Pressure continuation reuses only equilibria satisfying both GS and q
tolerances. The pressure-compensated initial FFprime prediction changes
the starting iterate, not the final convergence requirements.
When a pressure step is clipped to the target and fails, continuation now
halves the actual attempted pressure gap. This avoids repeating the same
failed target with an unchanged predictor.

For `I=0.925` and `I=1.0`, some points stalled with maximum q errors near
0.3% on the base mesh. Smaller finite-difference steps, a finer q-contour
grid alone, and alternative converged initial profiles did not recover
them. Increasing the native GS mesh to `NR=64, NZ=96`, together with
`contourNR=400`, recovered `(I,beta_0)=(0.925,0.06)` and `(1.0,0.05)`
with maximum constraint errors `0.02210%` and `0.02320%`, respectively.
Remaining points in these columns use that finer GS/q resolution.
The accepted physical pressure/q profiles and tolerances are unchanged.
Actual controls are saved in each recovered point's `solverConfiguration`.

For deep indentation, geometry recovery uses 4096 surface samples, 64
Fourier modes, and three projection iterations. This recovered the
`I=0.975` outer-surface closure failure without relaxing the 5% input
closure limit or removing `psi_N=0.95`. Full actual geometry controls are
saved in `point.ballooning.options`. The CSV includes these numerical
resolutions and the maximum measured input closure defect.

The GS source assembly reuses pressure/FFprime interpolants within each
assembly instead of reconstructing them in every element. Reference and
reused source vectors agreed exactly for both signed flux orientations,
row/column profiles, and quadrature degrees 7 and 10 on a native C-shape
mesh. The measured source-assembly speedup was 3.29x. The element weak form,
interpolation/extrapolation rules, local ordering, and quadrature are
unchanged. `verification/test_rhs_profile_reuse.m` repeats this check.

Observed eigenvalue sign changes are refined with additional complete
equilibrium/stability solves until the pressure bracket width is below
0.5% (at most 12 bisections). No second stable region is assumed.

For the 12-basis dense equilibria, increasing 96 to 192 cells/turn retained
the instability classification at `(I,beta_av)=(0.6,3.3851%),
(0.775,3.2876%),(0.6,13.2582%)`. At the first of these, maximum eigenvalues
were `2.15590e-6` and `2.18918e-6`. Detailed comparisons are saved in
`tmp/cshape_dense/dense_angular_resolution.mat`.

The map shows the current solver's sampled stability results. The preceding
Chance PRL benchmark did not establish agreement at high beta, and this
map is not an independent validation of the ballooning formulation. The
17-by-22 grid increases parameter sampling; it does not itself establish
GS mesh or radial/phase convergence.

## Run and resume

From `ballooning_equation` in MATLAB:

```matlab
addpath('examples');
run_dense_indented_cshape_beta_scan;
[scan,trace,summary] = collect_dense_indented_cshape_beta_scan;
```

For independent MATLAB processes, run disjoint columns:

```matlab
run_dense_indented_cshape_beta_scan(IndentationIndices=1:4:17);
run_dense_indented_cshape_beta_scan(IndentationIndices=2:4:17);
run_dense_indented_cshape_beta_scan(IndentationIndices=3:4:17);
run_dense_indented_cshape_beta_scan(IndentationIndices=4:4:17);
```

Use `-singleCompThread` for these independent CLI processes. Rerunning the
same command resumes completed point checkpoints without recalculating
them. A configuration mismatch raises an error; existing data are not
silently overwritten.

`SeedDirectory="..."` optionally loads converged native equilibria from
another set of column checkpoints as initial guesses. Every new point is
solved and tested against the current configuration; old stability values
are never imported. The 8-basis partial run is preserved in
`output/indented_cshape_dense_basis8_partial/`.

Default output: `ballooning_equation/output/indented_cshape_dense/`.
`column_XX.mat` contains each point, native equilibrium, predictor,
convergence diagnostics and full radial/phase eigenvalue array.
`trace_XX.mat` contains the refined transition brackets.
The collector writes `beta_scan.mat`, `samples.csv`, `marginal_trace.mat`,
`beta_av_indentation.png/.fig`, and `beta_av_construction_I.png`.
`critical_beta_indentation.png/.fig` and `critical_beta_construction_I.png`
show the first instability brackets alone, with their actual beta bracket
widths as bars. Markers in those boundary-only figures are bracket
midpoints, not additional equilibrium solves.
Only actual solves populate the table and plots; no interpolation fills
failed points.

After stopping any process writing the relevant column, the recovery
helpers can resume independently:

```matlab
out = fullfile(pwd,'output','indented_cshape_dense');
complete_dense_indented_cshape_column(fullfile(out,'column_14.mat'), ...
    NeighborFile=fullfile(out,'column_15.mat'),RefineMesh=true);
complete_dense_indented_cshape_column(fullfile(out,'column_17.mat'), ...
    NeighborFile=fullfile(out,'column_16.mat'),RefineMesh=true);
refine_dense_indented_cshape_geometry(fullfile(out,'column_16.mat'));
[scan,trace,summary] = collect_dense_indented_cshape_beta_scan;
```

Neighbor equilibria seed new solves only. Recovery preserves the base
configuration signature, while recording numerical overrides per point
and per refined trace. `marginal_brackets.csv` exposes the actual critical
beta brackets in addition to the full MAT trace.
