"""Select a cataloged bean equilibrium, inspect it and write a plain gfile copy."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from eq_process import Equilibrium, check_cocos
from equilibrium_io import write_geqdsk


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indentation", type=float, default=.25)
    parser.add_argument("--p0-kpa", type=float, default=925.)
    parser.add_argument("--dataset", type=Path, default=Path("data/bean"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.p0_kpa):
        parser.error("p0-kpa must be finite")
    with (args.dataset/"catalog.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    families = sorted({float(r["indentation"]) for r in rows})
    family = [r for r in rows if np.isclose(float(r["indentation"]), args.indentation, atol=1e-12, rtol=0)]
    if not family:
        parser.error(f"Choose an existing indentation: {families}")
    if args.output.exists():
        parser.error("Output already exists; choose a new directory")
    row = min(family, key=lambda r: abs(float(r["p0_kPa"])-args.p0_kpa))
    source = args.dataset/row["gfile"]
    meta = json.loads(Path(str(source)+".json").read_text())
    eq = Equilibrium.from_geqdsk(source)
    d = eq.data
    if check_cocos(d) != 1:
        raise RuntimeError("Expected COCOS 1")
    args.output.mkdir(parents=True)
    (args.output/"selected_catalog_row.json").write_text(json.dumps(row, indent=2)+"\n")
    # Rewriting through the writer updates the hash for the uncompressed bytes.
    # Preserve all original provenance and the original gzip transport hash.
    meta["source_compressed_gfile_sha256"] = meta["geqdsk_sha256"]
    plain = write_geqdsk(d, args.output/"selected.geqdsk", metadata=meta)
    fig, ax = plt.subplots(figsize=(6, 6))
    q_checks = []
    for level in (.1, .3, .5, .7, .9):
        surface = eq.extract_flux_surface(level, npoints=1024)
        integrals = eq.compute_contour_integrals(surface)
        q_checks.append(dict(psi_n=level, q_file=integrals.q_geqdsk,
                             q_from_contour=integrals.q_int, relative_error=integrals.q_relerr))
        ax.plot(surface.R, surface.Z, lw=1, label=f"psi_N={level:g}")
        if level == .5:
            field = eq.fields.eval(surface.R, surface.Z)
            # COCOS 1; EquilibriumFields itself also exposes psiR/psiZ directly.
            br, bz, bphi = field.psiZ/field.R, -field.psiR/field.R, field.Bphi
            np.savetxt(args.output/"surface_psiN_0p5.csv",
                np.column_stack((surface.R, surface.Z, br, bz, bphi, field.p)),
                delimiter=",", header="R_m,Z_m,BR_T,BZ_T,Bphi_T,p_Pa", comments="")
    ax.plot(np.r_[d.rbbbs, d.rbbbs[0]], np.r_[d.zbbbs, d.zbbbs[0]], "k-", lw=1.5, label="LCFS")
    ax.plot(d.rmaxis, d.zmaxis, "r+", label="magnetic axis")
    ax.set(xlabel="R [m]", ylabel="Z [m]", aspect="equal",
           title=f"{row['indent_id']}/{row['case']}; i={float(row['indentation']):g}")
    ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(args.output/"flux_surfaces.png", dpi=150); plt.close(fig)
    report = dict(selected=row, source=str(source), plain_gfile=str(plain),
        grid=[d.nw, d.nh], current_MA=d.current/1e6, q_axis=float(d.qpsi[0]),
        q_edge=float(d.qpsi[-1]), q_contour_checks=q_checks,
        note="Data inspection, not a new GS solve or whole-plasma stability test")
    (args.output/"inspection.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
    print("PASS: selected catalog equilibrium -> read/plot -> plain gfile")


if __name__ == "__main__":
    main()
