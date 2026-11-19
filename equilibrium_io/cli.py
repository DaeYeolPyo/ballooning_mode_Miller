"""Single DESC conversion and cataloged Chance-equilibrium batch export."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from .geqdsk import read_geqdsk, write_geqdsk


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _convert_case(job):
    from .desc import desc_to_data
    source, output, indentation, nr, include_native, resume, padding = job
    source, output = Path(source), Path(output)
    shape, case = source.parents[2].name, source.parent.name
    target = output/"gfiles"/shape/(case+".geqdsk.gz")
    native = output/"native"/shape/case/"equilibrium.h5"
    sidecar = Path(str(target)+".json")
    sha = _hash(source)
    if resume and target.is_file() and sidecar.is_file():
        meta = json.loads(sidecar.read_text())
        if (meta.get("source_h5_sha256") != sha or meta.get("grid", [None])[0] != nr
                or meta.get("geqdsk_sha256") != _hash(target)
                or meta.get("grid_padding_fraction", 0.) != padding):
            raise ValueError(f"Resume input/settings/hash mismatch: {target}")
        data = read_geqdsk(target)
    else:
        data, meta = desc_to_data(source, nr=nr, padding=padding)
        meta["source_relative_path"] = f"{shape}/equilibria/{case}/equilibrium.h5"
        meta["indentation"] = indentation
        write_geqdsk(data, target, metadata=meta)
        meta = json.loads(sidecar.read_text())
    if meta["F_surface_relative_std_max"] > 1e-3:
        raise ValueError(f"R*Bphi is insufficiently flux-constant: {source}")
    if include_native:
        native.parent.mkdir(parents=True, exist_ok=True)
        if native.exists() and _hash(native) != sha:
            raise ValueError(f"Different native equilibrium already exists: {native}")
        if not native.exists():
            shutil.copy2(source, native)
        shutil.copy2(source.with_name("metadata.json"), native.with_name("metadata.json"))
    original = meta["original_metadata"]
    return dict(indent_id=shape, indentation=indentation, case=case,
        p0_kPa=original["actual_p0_kPa"], beta_av_percent=original["beta_av_percent"],
        Ip_A=data.current, q_axis=float(data.qpsi[0]), q_edge=float(data.qpsi[-1]),
        F_edge_Tm=float(data.fpol[-1]), R_axis_m=data.rmaxis, Z_axis_m=data.zmaxis,
        psi_span_Wb_per_rad=data.sibry-data.simag, nr=data.nw, nz=data.nh,
        grid_padding_fraction=meta["grid_padding_fraction"],
        force_over_gradp_l2=original["force_over_gradp_l2"],
        F_surface_relative_std_max=meta["F_surface_relative_std_max"],
        inverse_map_max_residual_m=meta["inverse_map_max_residual_m"],
        Ip_relative_difference_from_DESC=meta.get("Ip_relative_difference_from_DESC"),
        accepted=original["accepted"], quality_passed=original["quality_passed"],
        gfile=str(target.relative_to(output)).replace("\\", "/"),
        native_h5=str(native.relative_to(output)).replace("\\", "/") if include_native else "",
        source_h5_sha256=sha, gfile_sha256=meta["geqdsk_sha256"])


def export_batch(source, output, *, nr=129, workers=1, include_native=False,
                 resume=False, padding=0.02):
    source, output = Path(source).resolve(), Path(output).resolve()
    if isinstance(workers, bool) or not isinstance(workers, int) or not 1 <= workers <= 32:
        raise ValueError("workers must be in 1..32")
    if not source.is_dir() or output == source or output.is_relative_to(source):
        raise ValueError("Output must be separate from the source equilibrium tree")
    paths = sorted(source.glob("indent_*/equilibria/p_*/equilibrium.h5"))
    if not paths:
        raise ValueError("No indent_*/equilibria/p_*/equilibrium.h5 found")
    jobs, plans = [], {}
    for path in paths:
        shape = path.parents[2].name
        run = path.parents[2]/"run.json"
        if shape not in plans:
            plans[shape] = json.loads(run.read_text())
        p = plans[shape]
        index = int(path.parent.name.split("_")[1])
        expected_pressure = p["pressure_kPa"][index]
        original = json.loads(path.with_name("metadata.json").read_text())
        if (not np.isclose(original["target_p0_kPa"], expected_pressure, rtol=0, atol=1e-9)
                or original["physics_id"] != p["physics_id"]):
            raise ValueError(f"Case does not match run physics: {path}")
        jobs.append((str(path), str(output), p["physics"]["indentation"], nr, include_native, resume, padding))
    for shape, plan in plans.items():
        actual = {p.parent.name for p in paths if p.parents[2].name == shape}
        expected = {f"p_{i:04d}" for i in range(len(plan["pressure_kPa"]))}
        if actual != expected:
            raise ValueError(f"Incomplete pressure set: {shape}; missing {sorted(expected-actual)}")
    output.mkdir(parents=True, exist_ok=True)
    provenance = output/"provenance"
    provenance.mkdir(exist_ok=True)
    for name in ("manifest.json", "stability_summary.csv", "surface_summary.csv", "summary.json"):
        if (source/name).is_file():
            shutil.copy2(source/name, provenance/("source_"+name))
    (provenance/"run_definitions.json").write_text(json.dumps(plans, indent=2)+"\n")
    rows, errors = [], []
    print(f"Exporting {len(jobs)} equilibria; COCOS 1, Wb/rad; exterior is not vacuum", flush=True)
    if workers == 1:
        for i, job in enumerate(jobs):
            rows.append(_convert_case(job))
            if (i+1) % 16 == 0:
                print(f"{i+1}/{len(jobs)}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_convert_case, job): job for job in jobs}
            for future in as_completed(futures):
                try:
                    rows.append(future.result())
                except Exception as error:
                    errors.append(dict(source=futures[future][0], error=str(error)))
                    if len(errors) <= 3:
                        print(f"Conversion failed: {errors[-1]}", flush=True)
                if (len(rows)+len(errors)) % 16 == 0:
                    print(f"{len(rows)} OK, {len(errors)} failed / {len(jobs)}", flush=True)
    if errors:
        (output/"conversion_errors.json").write_text(json.dumps(errors, indent=2)+"\n")
        raise RuntimeError(f"{len(errors)} conversions failed; see conversion_errors.json")
    rows.sort(key=lambda row: (row["indent_id"], row["case"]))
    with (output/"catalog.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    (output/"catalog.json").write_text(json.dumps(rows, indent=2, allow_nan=False)+"\n")
    (output/"conversion_errors.json").unlink(missing_ok=True)
    print(f"Completed {len(rows)} gfiles: {output}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    single = sub.add_parser("desc", help="Convert one supported DESC H5")
    single.add_argument("source", type=Path); single.add_argument("output", type=Path)
    single.add_argument("--nr", type=int, default=257)
    single.add_argument("--nz", type=int)
    single.add_argument("--padding", type=float, default=.02, help="Per-side fraction of LCFS R/Z extents")
    single.add_argument("--overwrite", action="store_true")
    batch = sub.add_parser("batch", help="Export a complete Chance run and catalog")
    batch.add_argument("source", type=Path); batch.add_argument("output", type=Path)
    batch.add_argument("--nr", type=int, default=129)
    batch.add_argument("--workers", type=int, default=1)
    batch.add_argument("--padding", type=float, default=.02)
    batch.add_argument("--include-native", action="store_true")
    batch.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.command == "desc":
        from .desc import export_desc_gfile
        print(export_desc_gfile(args.source, args.output, nr=args.nr, nz=args.nz,
                               padding=args.padding, overwrite=args.overwrite))
    else:
        export_batch(args.source, args.output, nr=args.nr, workers=args.workers,
                     include_native=args.include_native, resume=args.resume, padding=args.padding)


if __name__ == "__main__":
    main()
