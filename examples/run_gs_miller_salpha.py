"""One reference GS equilibrium -> fixed Miller geometry -> existing s-alpha scan.

Run from the project root: python -m examples.run_gs_miller_salpha
Requires the existing eq_process, local_equilibrium, salpha and
sturm_liouville modules on the Python import path.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import copy
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import time
from types import SimpleNamespace

import numpy as np
from threadpoolctl import threadpool_limits

from gs_fixed.solver import prepare_gs_problem, solve_gs, PicardOptions, QOptions, available_cpus
from gs_fixed.profiles import FFprimeModel, evaluate_ffprime_model
from gs_fixed.ballooning_adapter import build_ballooning_adapter
from salpha import prepare_fem_salpha_equilibrium, scan_salpha_diagram, plot_salpha_diagram


_WORKER_STATE = None
_WORKER_OPTIONS = None
_WORKER_LIMIT = None


def build_reference(args):
    """Edit this function to change the physical reference equilibrium."""
    r0, b0, a = 3., 3., .8
    kappa, delta = 1.5, .25
    theta = np.arange(120)*2*np.pi/120
    boundary = np.column_stack((
        r0+a*np.cos(theta+np.arcsin(delta)*np.sin(theta)),
        a*kappa*np.sin(theta),
    ))
    problem = prepare_gs_problem(
        boundary, r0=r0, b0=b0, a_minor=a, nr=args.nr, nz=args.nz,
        n_profile=101, q_options=QOptions(contour_nr=args.contour_nr),
        picard=PicardOptions(max_iterations=150, update_tolerance=1e-9,
                              residual_tolerance=1e-9),
    )
    u = problem.psi_n
    pprime = -1e6*(1-u**1.5)**2
    model = FFprimeModel(u, np.array([2., .2, -.1]), r0*b0)
    profile = evaluate_ffprime_model(model)
    reference = solve_gs(problem, pprime, ffprime=profile.ffprime)
    if not reference.converged:
        raise RuntimeError("Reference GS equilibrium did not converge.")
    definition = dict(r0=r0, b0=b0, a_minor=a, boundary_kappa=kappa,
                      boundary_delta=delta, ffprime_coefficients=model.coefficients.tolist(),
                      pprime_definition="dp/du=-1e6*(1-u**1.5)**2; u=psi_N", mode=1)
    return reference, definition


def _worker_init(state, options):
    global _WORKER_STATE, _WORKER_OPTIONS, _WORKER_LIMIT
    _WORKER_LIMIT = threadpool_limits(limits=1)
    _WORKER_STATE, _WORKER_OPTIONS = state, options


def _scan_block(rows, shat):
    result = scan_salpha_diagram(_WORKER_STATE, shat=shat,
                                 verbose=False, **_WORKER_OPTIONS)
    return rows, result


def run_scan(state, shat, alpha, theta0, args):
    options = dict(alpha=alpha, theta0=theta0,
                   theta_boundary=args.theta_max_pi*np.pi,
                   n_balloon_points=args.n_balloon_points,
                   enforce_periodicity=True, require_vector_consistency=True)
    # Each existing scan call requires at least two shat values.
    blocks = np.array_split(np.arange(len(shat)), max(1, len(shat)//2))
    cpus = available_cpus()
    workers = min(args.workers or cpus, cpus, len(blocks))
    print(f"s-alpha: {len(shat)} x {len(alpha)} x {len(theta0)}; workers={workers}", flush=True)
    started = time.perf_counter()
    if workers == 1:
        result = scan_salpha_diagram(state, shat=shat, verbose=True, **options)
    else:
        # Scans need only local geometry. Do not send the GS problem, its
        # SuperLU factor, Cartesian cache or P3 locator to scan workers.
        names = ("model", "targetPsiN", "localEquilibrium", "surface", "integrals",
                 "FFprime", "normalization", "shat", "alpha", "geometry")
        compact = SimpleNamespace(**{name:getattr(state, name) for name in names if hasattr(state, name)})
        result = None
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"),
                                  initializer=_worker_init, initargs=(compact, options)) as pool:
            futures = [pool.submit(_scan_block, rows, shat[rows]) for rows in blocks]
            for future in as_completed(futures):
                rows, part = future.result()
                if result is None:
                    result = copy(part)
                    result.shat = shat.copy()
                    # These arrays have a leading shat axis in the existing API.
                    names = ("lambdaMax", "lambdaTheta0", "bestTheta0", "failMessage",
                             "FFprimeUsed", "periodicityClosure", "validPhaseCount",
                             "qprime", "done", "rowSeconds")
                    for name in names:
                        a = getattr(part, name)
                        setattr(result, name, np.empty((len(shat),)+a.shape[1:], dtype=a.dtype))
                    # A block checkpoint signature is not a signature of the merged scan.
                    result.signature = None
                for name in names:
                    getattr(result, name)[rows] = getattr(part, name)
                print(f"  finished shat rows {rows[0]+1}-{rows[-1]+1}", flush=True)
        result.phaseComplete = result.validPhaseCount == len(theta0)
        result.numberFailed = int(np.count_nonzero(result.failMessage != ""))
        result.numberIncompletePoints = int(np.count_nonzero(~result.phaseComplete))
        result.scanElapsedSeconds = float(np.sum(result.rowSeconds))
        finite = result.lambdaMax[np.isfinite(result.lambdaMax)]
        result.lambdaRange = np.array([finite.min(), finite.max()]) if finite.size else np.full(2, np.nan)
        result.hasMarginalCrossing = bool(finite.size and finite.min() <= 0 <= finite.max())
    elapsed = time.perf_counter()-started
    return result, dict(requested_workers=args.workers, workers=workers,
                         available_cpus=cpus, blas_threads=1, seconds=elapsed)


def _save_reference(path, eq):
    np.savez_compressed(path, points=eq.mesh.points, elements=eq.mesh.elements,
        p1_points=eq.mesh.p1_points, p1_elements=eq.mesh.p1_elements,
        psi=eq.psi, psi_n=eq.psi_n, profile_psi_n=eq.problem.psi_n,
        pressure=eq.pressure, f=eq.f_profile.f, ffprime=eq.f_profile.ffprime,
        pprime=eq.pprime, q_psi_n=eq.q_result.psi_n, q=eq.q_result.q,
        axis_point=eq.axis_point, psi_axis=eq.state.psi_axis, dpsi=eq.state.dpsi)


def _save_miller(path, state):
    local, surface, p = state.localEquilibrium, state.surface, state.parameters
    np.savez_compressed(path, surface_r=surface.R, surface_z=surface.Z,
        theta=local.theta_geo, miller_r=local.R, miller_z=local.Z,
        bp=local.Bp, b2=local.B2, radial_fit_psi_n=p.radial_fit.psiN,
        radial_fit_r=p.radial_fit.r, radial_fit_R0=p.radial_fit.R0,
        radial_fit_delta=p.radial_fit.delta, radial_fit_kappa=p.radial_fit.kappa)


def _save_scan(path, result):
    complete = np.asarray(result.phaseComplete, dtype=bool)
    np.savez_compressed(path, shat=result.shat, alpha=result.alpha,
        theta0=result.theta0, theta=result.thetaGrid,
        lambda_theta0=result.lambdaTheta0,
        lambda_max=np.where(complete, result.lambdaMax, np.nan),
        lambda_max_observed=result.lambdaMax, best_theta0=result.bestTheta0,
        phase_complete=complete, valid_phase_count=result.validPhaseCount,
        fail_message=np.asarray(result.failMessage, dtype=str),
        ffprime_used=result.FFprimeUsed, periodicity_closure=result.periodicityClosure,
        qprime=result.qprime, pprime=result.pprime,
        equilibrium_shat=result.equilibriumShat, equilibrium_alpha=result.equilibriumAlpha)


def _source_records():
    import gs_fixed.solver, gs_fixed.ballooning_adapter
    import eq_process, local_equilibrium, salpha, sturm_liouville
    records = {}
    for module in (gs_fixed.solver, gs_fixed.ballooning_adapter,
                   eq_process, local_equilibrium, salpha, sturm_liouville):
        path = Path(module.__file__).resolve()
        records[module.__name__] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    path = Path(__file__).resolve()
    records["example"] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    return records


def render(state, result, args):
    if args.plot == "none":
        return
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6), constrained_layout=True)
    surface, local = state.surface, state.localEquilibrium
    ax.plot(np.r_[surface.R, surface.R[0]], np.r_[surface.Z, surface.Z[0]], label="GS surface")
    ax.plot(local.R, local.Z, "--", label="Miller fit")
    ax.set(xlabel="R [m]", ylabel="Z [m]", title=rf"Reference surface: $\psi_N={state.targetPsiN:.3f}$")
    ax.set_aspect("equal", adjustable="box")
    ax.legend()
    ax.grid(True, alpha=.25)
    figures = {"miller_fit":fig}
    if result is not None:
        plots = plot_salpha_diagram(result, show=False, mask_incomplete=True)
        figures.update(salpha=plots.stability, phase=plots.phase)
    if args.plot in ("save", "both"):
        directory = args.output/"figures"
        directory.mkdir(parents=True, exist_ok=True)
        for name, figure in figures.items():
            if figure is not None:
                target = directory/f"{name}.{args.plot_format}"
                figure.savefig(target, dpi=200, bbox_inches="tight")
                print(f"Saved figure: {target}", flush=True)
    if args.plot in ("show", "both"):
        plt.show(block=True)
    else:
        for figure in figures.values():
            if figure is not None:
                plt.close(figure)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nr", type=int, default=35)
    parser.add_argument("--nz", type=int, default=43)
    parser.add_argument("--contour-nr", type=int, default=200)
    parser.add_argument("--adapter-nr", type=int, default=257)
    parser.add_argument("--psi-n", type=float, default=.7)
    parser.add_argument("--surface-points", type=int, default=512)
    parser.add_argument("--fit-points", type=int, default=256)
    parser.add_argument("--volume-source", choices=("contour", "local"), default="contour")
    parser.add_argument("--shat-min", type=float, default=0.)
    parser.add_argument("--shat-max", type=float, default=5.)
    parser.add_argument("--n-shat", type=int, default=11)
    parser.add_argument("--alpha-min", type=float, default=0.)
    parser.add_argument("--alpha-max", type=float, default=6.)
    parser.add_argument("--n-alpha", type=int, default=13)
    parser.add_argument("--theta0-count", type=int, default=8)
    parser.add_argument("--theta-max-pi", type=float, default=5., help="theta boundary in units of pi")
    parser.add_argument("--n-balloon-points", type=int, default=151)
    parser.add_argument("--workers", type=int, default=1, help="scan processes; 0=allocation-aware auto")
    parser.add_argument("--prepare-only", action="store_true", help="save reference/Miller data without scanning")
    parser.add_argument("--plot", choices=("none", "show", "save", "both"), default="save")
    parser.add_argument("--plot-format", choices=("png", "pdf", "svg"), default="png")
    parser.add_argument("--output", type=Path, default=Path("results/gs_miller_salpha"))
    args = parser.parse_args()
    for name, minimum in (("nr", 5), ("nz", 5), ("contour_nr", 16), ("adapter_nr", 16),
                          ("surface_points", 32), ("fit_points", 16), ("n_shat", 2),
                          ("n_alpha", 2), ("theta0_count", 1), ("n_balloon_points", 21), ("workers", 0)):
        if getattr(args, name) < minimum:
            parser.error(f"--{name.replace('_','-')} must be >= {minimum}")
    if not 0 < args.psi_n < 1:
        parser.error("--psi-n must lie strictly between 0 and 1")
    if not np.isfinite(args.theta_max_pi) or args.theta_max_pi <= 0 or args.n_balloon_points % 2 == 0:
        parser.error("Require positive --theta-max-pi and odd --n-balloon-points")
    for name in ("shat", "alpha"):
        low, high = getattr(args, name+"_min"), getattr(args, name+"_max")
        if not np.all(np.isfinite([low, high])) or low >= high:
            parser.error(f"Require finite --{name}-min < --{name}-max")
    if args.plot != "none":
        import matplotlib
        if args.plot == "save":
            matplotlib.use("Agg", force=True)
        elif matplotlib.get_backend().lower() in ("agg", "pdf", "ps", "svg", "pgf", "cairo", "template"):
            parser.error("show/both needs a GUI backend; use --plot save on sbatch")
    return args


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    # Spawned workers inherit these before importing numerical libraries.
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    with threadpool_limits(limits=1):
        started = time.perf_counter()
        reference, definition = build_reference(args)
        gs_seconds = time.perf_counter()-started
        _save_reference(args.output/"reference_gs.npz", reference)
        print(f"Reference GS: {reference.mesh.n_dof} DOFs, {gs_seconds:.2f} s", flush=True)
        started = time.perf_counter()
        state = prepare_fem_salpha_equilibrium(
            reference, target_psiN=args.psi_n, adapter=build_ballooning_adapter,
            adapter_nr=args.adapter_nr, model="miller",
            n_surface_points=args.surface_points, n_fit_points=args.fit_points,
            radial_fit_half_width=.03, n_radial_fit=9, volume_source=args.volume_source,
        )
        miller_seconds = time.perf_counter()-started
        _save_miller(args.output/"miller_surface.npz", state)
        value = state.eqfunc.eval(state.surface.R, state.surface.Z)
        quality = dict(max_surface_flux_error=float(np.max(abs(value.psiN-args.psi_n))),
            raw_q_relative_error=float(abs(state.contourIntegrals.q_int/state.surface.q-1)),
            periodic_closure=float(state.thetaDiag.closure_used))
        print("Miller reference diagnostics:", json.dumps(quality), flush=True)
        p = state.parameters
        metadata = dict(python=platform.python_version(), numpy=np.__version__,
            arguments={key:str(v) if isinstance(v, Path) else v for key,v in vars(args).items()},
            definition=definition, sources=_source_records(),
            reference=dict(converged=reference.converged, n_dof=reference.mesh.n_dof,
                psi_n=state.targetPsiN, shat=state.shat, alpha=state.alpha,
                aN=state.normalization.aN, BN=state.normalization.BN,
                volume_source=state.geometry.volumeSource, V=state.integrals.V,
                Vprime=state.integrals.Vprime, q=state.surface.q,
                miller={name:float(getattr(p, name)) for name in (
                    "R0", "r", "delta", "kappa", "dR0_dr", "s_delta", "s_kappa")},
                quality=quality),
            timing=dict(gs_seconds=gs_seconds, miller_seconds=miller_seconds),
            scan=None,
        )
        result = None
        if not args.prepare_only:
            shat = np.linspace(args.shat_min, args.shat_max, args.n_shat)
            alpha = np.linspace(args.alpha_min, args.alpha_max, args.n_alpha)
            theta0 = 2*np.pi*(np.arange(args.theta0_count)-args.theta0_count//2)/args.theta0_count
            result, parallel = run_scan(state, shat, alpha, theta0, args)
            _save_scan(args.output/"salpha.npz", result)
            complete = result.numberIncompletePoints == 0
            metadata["scan"] = dict(complete=bool(complete), number_failed_phases=result.numberFailed,
                number_incomplete_points=result.numberIncompletePoints,
                n_shat=len(shat), n_alpha=len(alpha), n_theta0=len(theta0), parallel=parallel)
            print("Scan summary:", json.dumps(metadata["scan"]), flush=True)
        (args.output/"metadata.json").write_text(
            json.dumps(metadata, indent=2, allow_nan=False), encoding="utf-8")
        render(state, result, args)
    if result is not None and result.numberIncompletePoints:
        raise SystemExit("Incomplete scan: data and figures saved; inspect fail_message in salpha.npz")


if __name__ == "__main__":
    main()
