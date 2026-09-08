# Physics profile stage

This folder defines the flux normalization and source profiles used by the
fixed-boundary Grad-Shafranov equation.

## Convention

```text
psi_N = (psi - psi_axis)/(psi_boundary - psi_axis)
```

Thus the magnetic axis is `psi_N=0` and the LCFS is `psi_N=1`, regardless of
whether dimensional psi increases or decreases outward.

## Files

- `normalize_flux.m`: normalization with explicit extrapolate, clip, or error
  policy.
- `pprime.m`: constant, normalized-power, pressure-power, or tabulated
  derivative profile.
- `FFprime.m`: constant, normalized-power, F-power, or tabulated
  `F*dF/dpsi` profile.
- `evaluate_tabulated_profile.m`: monotone-grid interpolation on `psi_N`.
- `reconstruct_F.m`: reconstructs `F=R*Bphi`; a GEQDSK `FValues` table is
  used directly when supplied.
- `make_profiles.m`: creates and validates a consistent profile structure.
- `gs_source.m`: returns the FEM weak source
  `mu0*R*pprime + FFprime/R`.
- `verify_physics.m`: verifies analytic values and the connection to
  `fem/assemble_rhs.m`.

Run from `GS_fixed`:

```matlab
addpath('geometry', 'fem', 'physics')
report = verify_physics(true);
```

Linear verification profile:

```matlab
profiles = make_profiles(struct( ...
    'pprime',  struct('type','constant','value',Cp), ...
    'FFprime', struct('type','constant','value',CF)));
```

Normalized nonlinear derivative profile:

```matlab
profiles = make_profiles(struct( ...
    'flux', struct('psiAxis',psiAxis,'psiBoundary',0), ...
    'pprime', struct('type','normalized-power', ...
                     'axisValue',Cp,'exponent',2), ...
    'FFprime', struct('type','normalized-power', ...
                      'axisValue',CF,'exponent',1)));
```

GEQDSK-style tabulated profiles:

```matlab
profiles = make_profiles(struct( ...
    'flux',struct('psiAxis',simag,'psiBoundary',sibry), ...
    'pprime',struct('type','tabulated','psiN',psiN, ...
                    'values',pprimeValues,'method','pchip'), ...
    'FFprime',struct('type','tabulated','psiN',psiN, ...
                     'values',FFprimeValues,'method','pchip', ...
                     'FValues',FValues)));
```

During Picard iteration, `profiles.flux.psiAxis` should be updated from the
current iterate before evaluating a normalized profile. The upcoming solver
stage will own that update; this physics stage only evaluates the profiles.
