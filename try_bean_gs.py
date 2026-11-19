"""Colleague acceptance example: custom bean boundary -> GS -> gfile -> readback.

Run AFTER installing the repository, from its top-level directory.
The boundary follows Chance et al. (1983), Eq. (2), with requested d/(2a).
This example prescribes pressure and F*dF/dpsi_N; q is an output.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from threadpoolctl import threadpool_limits

from eq_process import Equilibrium, check_cocos
from gs_fixed.solver import prepare_gs_problem, solve_gs, PicardOptions, QOptions
from gs_fixed.post import build_post_data, compute_global_parameters, plot_gs_equilibrium


def periodic_minimum(function):
    t = np.linspace(0., np.pi, 257)
    values = function(t)
    ids = np.flatnonzero((values[1:-1] <= values[:-2])
                         & (values[1:-1] <= values[2:]))+1
    candidates = [float(values[0]), float(values[-1])]
    for i in ids:
        result = minimize_scalar(function, bounds=(t[i-1], t[i+1]),
                                 method="bounded", options={"xatol": 1e-13})
        if not result.success:
            raise RuntimeError("Boundary extremum search failed")
        candidates.append(float(result.fun))
    return min(candidates)


def bean_boundary(indentation, r0, a, kappa, n):
    """Indentation i=(R(pi)-R_min)/(R_max-R_min), not Miller delta."""
    x = lambda t, c: (1+.6*np.cos(t))*np.cos(c*np.sin(t))
    y = lambda t, c: (1+.6*np.cos(t))*np.sin(c*np.sin(t))
    def depth(c):
        xmin = periodic_minimum(lambda t: x(t, c))
        return (.4-xmin)/(1.6-xmin)
    cmax = np.pi-1e-5
    if not np.isfinite(indentation) or not 0 < indentation < depth(cmax):
        raise ValueError("Indentation is outside the positive, simple bean branch")
    c = brentq(lambda c: depth(c)-indentation, 0., cmax, xtol=1e-13)
    xmin, xmax = periodic_minimum(lambda t: x(t, c)), 1.6
    ymax = -periodic_minimum(lambda t: -y(t, c))
    raw_a, center = (xmax-xmin)/2, (xmax+xmin)/2
    t = np.linspace(0., 2*np.pi, n, endpoint=False)
    points = np.column_stack((r0+a*(x(t, c)-center)/raw_a,
                              kappa*a*y(t, c)/ymax))
    if np.any(points[:, 0] <= 0):
        raise ValueError("The boundary reaches nonpositive R")
    return points, dict(C=c, B=.6, measured_indentation=depth(c),
                        R_min_m=r0-a, R_max_m=r0+a, Z_max_m=kappa*a)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indentation", type=float, default=.25)
    parser.add_argument("--r0", type=float, default=3.2, help="Geometric major radius [m]")
    parser.add_argument("--a", type=float, default=.8, help="Geometric minor radius [m]")
    parser.add_argument("--kappa", type=float, default=1.386, help="Z_max/a")
    parser.add_argument("--b0", type=float, default=3., help="F_boundary/r0 [T]")
    parser.add_argument("--p0-kpa", type=float, default=50.)
    parser.add_argument("--ff-scale", type=float, default=-2., help="F*dF/dpsi_N = scale*(1-u) [T^2 m^2]")
    parser.add_argument("--mesh-nr", type=int, default=31)
    parser.add_argument("--mesh-nz", type=int, default=43)
    parser.add_argument("--contour-nr", type=int, default=256)
    parser.add_argument("--export-nr", type=int, default=257)
    parser.add_argument("--boundary-points", type=int, default=192)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    args = parser.parse_args()
    if (not np.all(np.isfinite([args.r0, args.a, args.kappa, args.b0, args.p0_kpa, args.ff_scale]))
            or min(args.r0-args.a, args.a, args.kappa, args.b0, args.p0_kpa) <= 0):
        parser.error("Require r0>a>0, kappa>0, b0>0, p0>0 and finite input")
    if min(args.mesh_nr, args.mesh_nz, args.boundary_points) < 16:
        parser.error("Mesh candidate grid and boundary counts must be >=16")
    if args.output.exists():
        parser.error("Output already exists; choose a new --output directory")
    boundary, shape = bean_boundary(args.indentation, args.r0, args.a,
                                    args.kappa, args.boundary_points)
    with threadpool_limits(limits=1):
        problem = prepare_gs_problem(boundary, r0=args.r0, b0=args.b0, a_minor=args.a,
            nr=args.mesh_nr, nz=args.mesh_nz, n_profile=201,
            picard=PicardOptions(max_iterations=250, omega=.5,
                update_tolerance=1e-9, residual_tolerance=1e-9, error_on_failure=True),
            q_options=QOptions(contour_nr=args.contour_nr))
        u = problem.psi_n
        # u = outward normalized poloidal flux. Inputs are d/du, not d/dpsi.
        pprime = -4*(args.p0_kpa*1000)*u*(1-u*u)
        ffprime = args.ff_scale*(1-u)
        eq = solve_gs(problem, pprime, ffprime=ffprime, p_edge=0.)
        if not eq.converged:
            raise RuntimeError("GS did not converge")
        global_data = compute_global_parameters(build_post_data(eq))
        args.output.mkdir(parents=True)
        np.savetxt(args.output/"boundary.csv", boundary, delimiter=",", header="R_m,Z_m", comments="")
        gfile = eq.export_gfile(args.output/"equilibrium.geqdsk", nr=args.export_nr)
        loaded = Equilibrium.from_geqdsk(gfile)
        cocos = check_cocos(loaded.data)
        if cocos != 1 or not np.all(np.isfinite(loaded.data.psirz)):
            raise RuntimeError("Invalid exported G-EQDSK data")
        surface_checks = []
        for level in (.1, .5, .9):
            surf = loaded.extract_flux_surface(level, npoints=2048)
            integrals = loaded.compute_contour_integrals(surf)
            surface_checks.append(dict(psi_n=level, q_file=integrals.q_geqdsk,
                q_from_contour=integrals.q_int, relative_error=integrals.q_relerr))
        # This is an interchange acceptance limit, not a physical-accuracy claim.
        max_q_error = max(abs(row["relative_error"]) for row in surface_checks)
        passed = bool(max_q_error < .01)
        report = dict(passed=passed, converged=bool(eq.converged), cocos=cocos,
            requested=dict(indentation=args.indentation, r0_m=args.r0, a_m=args.a,
                kappa=args.kappa, b0_T=args.b0, p0_kPa=args.p0_kpa, ff_scale=args.ff_scale),
            boundary=shape, pressure_axis_actual_kPa=float(eq.pressure[0]/1000),
            q_axis=float(eq.q_result.q[0]), q_edge=float(eq.q_result.q[-1]),
            iterations=eq.inner_result.iterations,
            final_update=float(eq.inner_result.update_history[-1]),
            final_residual=float(eq.inner_result.residual_history[-1]),
            q_readback_checks=surface_checks, q_readback_max_relative_error=max_q_error,
            global_parameters=global_data,
            gfile="equilibrium.geqdsk", exterior="normal continuation; not a vacuum solution")
        (args.output/"acceptance.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        fig, _ = plot_gs_equilibrium(eq)
        fig.savefig(args.output/"equilibrium.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
    print(json.dumps(report, indent=2))
    if not passed:
        raise RuntimeError("Gfile q readback differs by >=1%; refine mesh/contour/export grids")
    print("PASS: new bean GS -> gfile -> readback")


if __name__ == "__main__":
    main()
