# Gfile export definitions

The byte layout follows the [GA G-EQDSK specification](https://fusion.gat.com/conferences/snowmass/working/mfe/physics/p3/equilibria/g_eqdsk_s.pdf):
48-character label, three I4 integers, E16.9 values in groups of five, and an
R-fast Fortran-order rectangular flux array. Profiles use `nw` samples uniform
in normalized **poloidal** flux. R,Z are metres, p is Pa, F is T m, Ip is A.

## Coordinate convention

Adapters export [COCOS 1](https://www.epfl.ch/research/domains/swiss-plasma-center/wp-content/uploads/2018/10/Sauter_COCOS_Tokamak_Coordinate_Conventions.pdf):

```text
psi(axis)=0; psi(boundary)=positive span; psi units=Wb/rad
BR = +psi_Z/R; BZ = -psi_R/R; Bphi = F/R
q has the sign of F for this canonical positive-current orientation
```

Gfiles cannot encode the complete coordinate convention themselves. Keep the
adjacent JSON sidecar with the file. The generic writer preserves the supplied
signs; it does not perform a COCOS conversion for arbitrary input data.

## GS export

```python
eq.export_gfile("equilibrium.geqdsk", nr=257, nz=385)
# or: from gs_fixed.geqdsk import export_gfile
```

Let `d = state.dpsi`, `a = state.psi_axis`:

```text
psi_out = sign(d)*(psi_native-a)
pprime_out = (dp/dpsi_N)/abs(d)
ffprime_out = (F dF/dpsi_N)/abs(d)
Ip_out = -sign(d)*Ip_native
Bp_out = -sign(d)*Bp_native
```

The nodal-axis normalization used by the solver is preserved. Its difference
from the recovered continuous P3 axis flux is reported as
`axis_flux_offset_Wb_per_rad`; the exporter does not silently change the
profiles' normalization. Interior flux values and gradients come from the
original P3 triangles. Source derivatives, p and F preserve the solver's linear
profile interpolation. q uses the existing postprocessor interpolation, with
its extrapolated edge value explicitly identified in metadata.

## DESC H5 export

```bash
python -m equilibrium_io desc equilibrium.h5 output.geqdsk --nr 257
python -m equilibrium_io batch CHANCE_RUN data/bean --nr 129 --workers 4 --include-native
```

For DESC toroidal flux `Phi` and toroidal-flux radius `rho`:

```text
dpsi/drho = Phi*rho*iota(rho)/pi
psi(rho)  = Phi/pi * integral_0^rho t*iota(t) dt
q = 1/iota
pprime = p_rho / psi_rho
F = <R Bphi>_surface, weighted by DESC sqrt(g)
FFprime = F*F_rho/psi_rho
```

The rho*iota integral is the exact piecewise-polynomial integral of the saved
Hermite spline. Uniform poloidal-flux samples are inverted with a bracketing
root solver. F and its radial derivative use analytic Zernike derivatives;
the axis limits use even-in-rho extrapolation. Current is independently
integrated from the boundary poloidal field via Ampere's law. The small finite
angular variation of R Bphi and current difference from DESC metadata are
reported. The reference B is `F_edge/R_geo`, rather than the input vacuum B
being copied over a finite-beta equilibrium with a different edge F.

Rectangular grid points are inverted to rho,theta using Newton iteration on
Fourier geometry. A radial amplitude spline accelerates the inversion; its
sampled error is checked against independent exact Zernike evaluation. Mapping
residuals are recorded. This is not scattered interpolation of psi labels.
Only the archive representation validated for DESC 0.17.3 is supported.

## Exterior and accuracy

Fixed-boundary solutions do not contain a physical exterior magnetic field.
To provide a finite rectangular gfile grid, outside points are filled with
`psi_b + |grad(psi)_boundary| * distance_to_boundary`, using the nearest LCFS
segments and interpolated boundary gradient magnitudes. This continuation
matches the first normal derivative locally but does **not** solve the vacuum
GS equation. There are no coils or limiter/wall data in these files.

Use the supplied LCFS to distinguish the supported interior from the exterior.
Vacuum/external-mode calculations require another exterior solution. The grid
also interpolates the original spectral/P3 interior; higher derivatives and
ballooning growth rates require their own resolution checks. The stored
spectral stability labels are not certification of a gfile round-trip.

`.gz` changes transport only. `gzip -dk file.geqdsk.gz` produces a conventional
plain-text gfile. Both repository readers accept plain and gzip files. Existing
destinations are refused unless `overwrite=True` / `--overwrite` is explicit.
