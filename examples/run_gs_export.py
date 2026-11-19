"""Solve one shaped fixed-boundary equilibrium and write a gfile + metadata."""
import argparse
from pathlib import Path
import numpy as np

from gs_fixed.solver import prepare_gs_problem, solve_gs, PicardOptions, QOptions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/gs_example"))
    parser.add_argument("--nr", type=int, default=257, help="R grid size in exported gfile")
    args = parser.parse_args()
    theta = np.linspace(0, 2*np.pi, 120, endpoint=False)
    boundary = np.column_stack((3+.8*np.cos(theta+np.arcsin(.25)*np.sin(theta)),
                                1.2*np.sin(theta)))
    problem = prepare_gs_problem(boundary, r0=3., b0=3., a_minor=.8,
        nr=25, nz=31, n_profile=201,
        picard=PicardOptions(max_iterations=150), q_options=QOptions(contour_nr=240))
    u = problem.psi_n
    eq = solve_gs(problem, -1e6*(1-u**1.5)**2, ffprime=2*(1-u))
    path = eq.export_gfile(args.output/"equilibrium.geqdsk", nr=args.nr)
    print(f"converged={eq.converged}; {path}")


if __name__ == "__main__":
    main()
