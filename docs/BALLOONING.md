# Ballooning solver routes

1. `sturm_liouville.Sturm_Liouville`: reusable P2 finite-element solver for
   `(g X')'+c X=lambda f X`, with its existing boundary conditions and largest
   algebraic eigenvalue selection. Coefficient arrays and theta grid are inputs.
2. `salpha.py`: local-equilibrium preparation, s-alpha scans, marginal roots,
   and plots. Models include direct flux-surface geometry, Miller, circular,
   C-shape and indented Miller. The GS `ballooning_adapter` reads original FEM
   values/gradients; gfile input uses the rectangular-grid reader.
3. `ballooning.desc`: the existing Chance-study spectral H5 coefficient path
   and finite-difference eigenproblem in `zeta=q*theta_PEST`. It retains the
   original geometric-curvature arithmetic and Dirichlet endpoints. The
   force-balance curvature route is available as a diagnostic coefficient but
   is not silently substituted into `solve_desc_ballooning`.

Positive largest lambda denotes instability in all routes. Numerical
normalizations differ between routes: compare dimensional definitions and
normalization before comparing magnitudes. A largest eigenvalue for one
surface/phase is not a global plasma-stability verdict.

```python
import numpy as np
from sturm_liouville import Sturm_Liouville

t = np.linspace(-np.pi/2, np.pi/2, 65)
result = Sturm_Liouville(t, np.ones_like(t), np.zeros_like(t),
                        np.ones_like(t), numEigenvalues=1)
print(result.eigval[0])  # approaches -1, the analytic Dirichlet result
```

For bean cases, the original H5 is the preferred source for derivative-sensitive
checks. The catalog's 129-point R gfiles are an interchange dataset. Re-export
at larger resolution, increase domain/cells and scan additional flux surfaces
and phases for a new stability study. The original Chance scan used finite
samples; its metadata records what was actually checked.
