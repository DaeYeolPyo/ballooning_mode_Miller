# Post-processing stage

This folder diagnoses the P1 fixed-boundary equilibrium without interpolating
through the empty notch of a concave C-shaped plasma.

## Files

- `find_axis.m`: interior flux extremum plus an accepted/rejected local
  quadratic refinement of `(R_axis,Z_axis)`.
- `extract_flux_contours.m`: marching-triangle level-set extraction and
  segment stitching directly on the plasma mesh.
- `find_critical_points.m`: discrete P1 O- and X-point candidates from cyclic
  neighbor signs.
- `check_topology.m`: one closed axis-enclosing contour per normalized level,
  pairwise nesting, and optional single-O/no-X requirements.
- `compute_B.m`: element-exact and recovered nodal `B_R`, `B_Z`, with optional
  `B_phi=F/R`.
- `compute_magnetic_curvature.m`: evaluates the full cylindrical
  `kappa=(b dot grad)b`, `grad(p)=pprime(psi) grad(psi)`, their signed dot
  product, and area/volume/normalized-flux-band statistics.
- `plot_curvature.m`: a signed, zero-centered map of `kappa dot grad(p)`
  with every zero-level boundary overlaid as a black dashed curve.
- `plot_flux.m`: mesh-safe colored flux plot, contours, LCFS, and axis marker.
- `verify_post.m`: verifies all diagnostics on a concave C-shaped equilibrium
  and checks `compute_B` against the manufactured field from `psi=R`.

Run from `GS_fixed`:

```matlab
addpath('geometry','fem','physics','solver','post')
report = verify_post(true);
```

Typical use:

```matlab
axisData = find_axis(mesh,psi,struct('mode','min'));
topology = check_topology(mesh,psi,axisData, ...
    struct('failOnFailure',true));
field = compute_B(mesh,psi,F_of_psi);
curvature = compute_magnetic_curvature(mesh,psi,field,profiles, ...
    struct('psiAxis',axisData.psi,'boundaryValue',0));
plot_flux(mesh,psi,struct('axisData',axisData));
plot_curvature(mesh,curvature);
```

The curvature formula includes the cylindrical-basis terms `-b_phi^2/R`
and `b_R*b_phi/R`. Element values are the primary diagnostic because
curvature differentiates a recovered P1 magnetic field; nodal values are
area-weighted visualization recoveries. The code reports the signed value
of `kappa dot grad(p)` without relabeling either sign as good or bad. Its
physical interpretation depends on the pressure and flux orientation used
by the equilibrium.

## Cosine-boundary curvature scan

Run the nonlinear screening scan from the project root with:

```matlab
result = scan_cosine_curvature();
result.ranking(1:10,:)
```

The scan sets `kappa` separately for every candidate so the actual
bounding-box elongation remains fixed, then rejects invalid meshes,
off-midplane axes, multiple O-points, X-points, and failed nested-contour
checks. With the `main_cos` flux/profile convention, positive
`kappa dot grad(p)` is counted as bad curvature. The ranking uses the
positive fraction of sign-active toroidal volume; near-zero regions are
reported separately so a candidate cannot appear improved merely because
the pressure gradient vanishes over more of the edge.

The contour test is the primary topology result. `find_critical_points` is a
piecewise-linear vertex diagnostic and should be interpreted together with
mesh refinement. A refined continuous critical-point search can be added when
the solver is upgraded from P1 to P2.
