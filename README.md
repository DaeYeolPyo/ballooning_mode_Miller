# Grad-Shafranov equilibrium and ballooning solvers

Python fixed-boundary P3 finite-element Grad-Shafranov solver, ideal ballooning
solvers, and a catalog of 344 bean-shape equilibria from the Chance Fig. 3 study.
The earlier MATLAB implementation remains available in Git history.

## Install

Python 3.10 or newer:

```bash
python -m pip install -e .
```

The core uses NumPy, SciPy, Matplotlib, ContourPy, Triangle, threadpoolctl and
h5py. Reading the included DESC 0.17.3 archives and exporting gfiles does not
require DESC or JAX. The H5 adapter supports axisymmetric Fourier-Zernike
geometry with HermiteSplineProfile pressure/iota; other inputs are rejected.

## Solve GS and export a gfile

```bash
python -m examples.run_gs_export --output results/gs_example
```

```python
from gs_fixed.solver import prepare_gs_problem, solve_gs

# boundary_points: ordered (N,2) array of R,Z coordinates in metres
problem = prepare_gs_problem(boundary_points, r0=3., b0=3., a_minor=.8)
u = problem.psi_n
eq = solve_gs(problem, -1e6*(1-u**1.5)**2, ffprime=2*(1-u))
eq.export_gfile("equilibrium.geqdsk", nr=257)
```

`pprime` and `ffprime` solver inputs mean **dp/dpsi_N** and **F dF/dpsi_N**.
Export converts them to dimensional derivatives per Wb/rad. Only converged
equilibria are exported; q is computed if it was skipped during the solve.
The prescribed-q path and serial/process round-trip example are retained:

```bash
python -m examples.run_gs_roundtrip --workers 1
```

## Ballooning

`sturm_liouville.Sturm_Liouville` solves the P2 finite-element generalized
eigenproblem `(g X')' + c X = lambda f X`. Positive largest lambda denotes
instability. `salpha.py` prepares local geometry and scans s-alpha diagrams;
the GS adapter evaluates the original FEM field directly.

```bash
python -m examples.run_gs_miller_salpha --help
python -m examples.ballooning_from_gfile data/bean/gfiles/indent_00/p_0000.geqdsk.gz --psi-n .5
```

The DESC archive path retains the original spectral geometry/curvature solver:

```python
from ballooning.desc import solve_desc_ballooning

result = solve_desc_ballooning(
    "data/bean/native/indent_00/p_0000/equilibrium.h5", .5,
    theta0=0., theta_max_pi=16, cells_per_turn=512)
print(result["lambda"])
```

One nonpositive sampled eigenvalue does not establish global stability. Check
surface/phase sampling and domain/grid convergence before interpreting a scan.
See [ballooning guide](docs/BALLOONING.md) and [s-alpha API](docs/salpha_usage.md).

## Bean equilibria for follow-up work

[Dataset guide](data/bean/README.md) and [catalog CSV](data/bean/catalog.csv)
identify the eight indentation families, pressure points, beta, current,
quality metrics, gfile paths, native H5 paths and SHA-256 hashes.

```python
from eq_process import Equilibrium

eq = Equilibrium.from_geqdsk("data/bean/gfiles/indent_00/p_0000.geqdsk.gz")
surface = eq.extract_flux_surface(.5)
```

To produce a higher-resolution file from an included native equilibrium:

```bash
python -m equilibrium_io desc data/bean/native/indent_00/p_0000/equilibrium.h5 refined.geqdsk --nr 513
```

All new gfiles use **COCOS 1**, reduced poloidal flux in **Wb/rad**, pressure
in Pa, F in T m, and current in A. Each has a mandatory provenance `.json`.
The sampling box includes a 2% margin on each side of the LCFS bounding box,
so midplane readers can reach the boundary flux before hitting the grid edge.
The supplied equilibria are fixed-boundary: **the exterior grid is a normal
linear continuation, not a vacuum solution**. Use the files inside the supplied
LCFS. External-mode, coil, wall, or vacuum-region calculations need a physical
exterior solution. The limiter contour is empty; the LCFS is not relabeled as a wall.

The distributed 129-point R grids support compact interchange. Re-export and
check convergence for derivative-sensitive work, especially ballooning near
the edge. Native H5 evaluation avoids this rectangular-grid interpolation.

## Validation

```bash
python -m pip install -e ".[test]"
python -m unittest test.test_geqdsk test.test_desc_export test.test_ballooning
python -m test.test_fem
python -m test.test_mesh
python -m test.test_assembly
python -m test.test_profiles
python -m test.test_post
python -m test.test_post_global
python -m test.test_solver
python scripts/validate_bean_catalog.py data/bean
```

See [export definitions](docs/GFILE_EXPORT.md), [source provenance](docs/SOURCE_IMPORT.md)
and [dataset validation](data/bean/validation.json). Examples write to ignored
`results/`; the packaged wheel contains code, while the Git repository contains
the equilibrium dataset.
