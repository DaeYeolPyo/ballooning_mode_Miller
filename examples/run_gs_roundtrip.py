"""Mode-1 -> prescribed-q round trip and serial/process benchmark."""
import argparse
import json
import os
import platform
from pathlib import Path
import time
import numpy as np
from threadpoolctl import threadpool_limits
from gs_fixed.solver import (
    prepare_gs_problem,solve_gs,PicardOptions,PQOptions,ParallelOptions,
    QOptions,available_cpus,
)
from gs_fixed.fem import assemble_rhs
from gs_fixed.profiles import FFprimeModel,QTarget,evaluate_ffprime_model
from gs_fixed.post import (
    plot_gs_equilibrium,
    build_post_data,
    compute_global_parameters,
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers",type=int,default=0,help="0=allocation-aware automatic; 1=serial")
    parser.add_argument("--nr",type=int,default=35)
    parser.add_argument("--nz",type=int,default=43)
    parser.add_argument("--contour-nr",type=int,default=160)
    parser.add_argument("--basis",type=int,default=3)
    parser.add_argument("--boundary-points",type=int,default=120)
    parser.add_argument("--benchmark",action="store_true",help="Run serial and process variants with startup included")
    parser.add_argument("--benchmark-rhs",action="store_true")
    parser.add_argument("--output",type=Path,default=Path("results/gs_roundtrip"))
    parser.add_argument(
        "--plot",
        choices=("none", "show", "save", "both"),
        default="none",
    )
    parser.add_argument("--plot-nr", type=int, default=240)
    parser.add_argument(
        "--plot-format", choices=("png", "pdf", "svg"), default="png",
    )

    args=parser.parse_args()

    if args.plot_nr < 16:
        parser.error("--plot-nr must be >=16")
    if args.plot != "none":
        import matplotlib
        if args.plot == "save":
            matplotlib.use("Agg", force=True)
        elif matplotlib.get_backend().lower() in (
            "agg", "pdf", "ps", "svg", "pgf", "cairo", "template",
        ):
            parser.error(
                "--plot show/both requires a GUI backend. "
                "Use --plot save on sbatch."
            )
    if not 1<=args.basis<=12: parser.error("--basis must be in 1..12")
    args.output.mkdir(parents=True,exist_ok=True)
    theta=np.arange(args.boundary_points)*2*np.pi/args.boundary_points
    boundary=np.column_stack((3+.8*np.cos(theta+np.arcsin(.25)*np.sin(theta)),
                               1.2*np.sin(theta)))
    with threadpool_limits(limits=1):
        problem=prepare_gs_problem(boundary,r0=3.,b0=3.,a_minor=.8,
            nr=args.nr,nz=args.nz,n_profile=101,
            q_options=QOptions(contour_nr=args.contour_nr),
            picard=PicardOptions(max_iterations=150,update_tolerance=1e-9,residual_tolerance=1e-9))
        s=problem.psi_n
        pp=-1e6*(1-s**1.5)**2
        exact_coefficients=np.zeros(args.basis); exact_coefficients[0]=2.
        if args.basis>1: exact_coefficients[1]=.2
        if args.basis>2: exact_coefficients[2]=-.1
        exact_profile=evaluate_ffprime_model(FFprimeModel(s,exact_coefficients,9.))
        start=time.perf_counter()
        reference=solve_gs(problem,pp,ffprime=exact_profile.ffprime)
        if not reference.converged: raise RuntimeError("Reference equilibrium failed.")
        reference_seconds=time.perf_counter()-start
        target=QTarget(reference.q_result.psi_n,reference.q_result.q)
        guess=np.zeros(args.basis); guess[0]=1.
        initial=FFprimeModel(s,guess,9.)
        options=PQOptions(max_iterations=15,q_relative_tolerance=1e-4,
            q_rms_log_tolerance=1e-4,regularization=0.,verbose=True)
        report=dict(python=platform.python_version(),numpy=np.__version__,
            cpus_available=available_cpus(),slurm_job_id=os.environ.get("SLURM_JOB_ID"),
            n_dof=problem.mesh.n_dof,n_elements=problem.mesh.n_elements,
            n_basis=args.basis,contour_nr=args.contour_nr,
            preparation_seconds=problem.preparation_seconds,reference_seconds=reference_seconds)
        if args.benchmark_rhs:
            start=time.perf_counter()
            for _ in range(10):
                cached=problem.assemble_rhs(problem.initial_psi,pp,exact_profile.ffprime,-.5)
            cached_seconds=(time.perf_counter()-start)/10
            start=time.perf_counter()
            uncached=assemble_rhs(problem.mesh,problem.initial_psi,problem.quad,pp,exact_profile.ffprime,-.5)
            loop_seconds=time.perf_counter()-start
            np.testing.assert_allclose(cached,uncached,rtol=5e-14,atol=2e-15)
            report["rhs"]=dict(cached_seconds=cached_seconds,element_loop_seconds=loop_seconds,
                               speedup=loop_seconds/cached_seconds)
        requested=[("serial",1),("parallel",args.workers)] if args.benchmark else [("run",args.workers)]
        runs={}
        for name,workers in requested:
            print(f"\n{name}: requested workers={workers}, DOFs={problem.mesh.n_dof}",flush=True)
            start=time.perf_counter()
            eq=solve_gs(problem,pp,q_target=target,ffprime_model=initial,outer_options=options,
                        parallel=ParallelOptions(max_workers=workers))
            elapsed=time.perf_counter()-start
            result=eq.pq_result
            report[name]=dict(seconds=elapsed,converged=eq.converged,
                termination_reason=result.termination_reason,outer_iterations=result.outer_iterations,
                function_evaluations=result.function_evaluations,
                max_relative_q=result.final_metrics.max_abs_relative,
                rms_log_q=result.final_metrics.rms_log,
                coefficients=result.final_model.coefficients.tolist(),
                parallel=result.parallel_info)
            np.savez_compressed(args.output/(name+".npz"),points=eq.mesh.points,
                elements=eq.mesh.elements,psi=eq.psi,psi_n=eq.psi_n,
                profile_psi_n=s,pressure=eq.pressure,f=eq.f_profile.f,g=eq.f_profile.g,
                ffprime=eq.f_profile.ffprime,q_psi_n=eq.q_result.psi_n,q=eq.q_result.q,
                coefficients=result.final_model.coefficients)
            runs[name]=eq
            if not eq.converged:
                (args.output/"summary.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
                raise RuntimeError(f"{name} failed: {result.termination_reason}")
        if args.benchmark:
            a,b=runs["serial"],runs["parallel"]
            np.testing.assert_allclose(a.psi,b.psi,rtol=2e-10,atol=2e-11)
            np.testing.assert_allclose(a.pq_result.final_model.coefficients,
                                       b.pq_result.final_model.coefficients,rtol=2e-9,atol=2e-9)
            report["overall_process_speedup"]=report["serial"]["seconds"]/report["parallel"]["seconds"]
            report["serial_parallel_max_psi_difference"]=float(np.max(abs(a.psi-b.psi)))
        (args.output/"summary.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
        print(json.dumps(report,indent=2),flush=True)

        if args.plot != "none":
            import matplotlib.pyplot as plt

            save = args.plot in ("save", "both")
            figure_dir = args.output / "figures"

            if save:
                figure_dir.mkdir(parents=True, exist_ok=True)

            for name, eq in runs.items():
                quantities = compute_global_parameters(build_post_data(eq))

                fig, _ = plot_gs_equilibrium(
                    eq,
                    nr=args.plot_nr,
                    q_target=target,
                    quantities=quantities,
                )

                if save:
                    path = figure_dir / f"{name}_equilibrium.{args.plot_format}"
                    fig.savefig(path, dpi=200, bbox_inches="tight")
                    print(f"Saved figure: {path}", flush=True)

                if args.plot == "save":
                    plt.close(fig)

            if args.plot in ("show", "both"):
                plt.show(block=True)


if __name__=="__main__":
    main()
