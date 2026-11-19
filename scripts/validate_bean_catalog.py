"""Verify every catalog file/hash and optionally parse all gfiles with FreeQDSK."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator, interp1d

# Also works as a script before the project has been installed.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from eq_process import read_geqdsk, check_cocos


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset", type=Path)
    parser.add_argument("--require-freeqdsk", action="store_true")
    args = parser.parse_args()
    try:
        from freeqdsk import geqdsk
    except ImportError:
        geqdsk = None
    if args.require_freeqdsk and geqdsk is None:
        raise RuntimeError("Install the project's [test] extra to use FreeQDSK")
    root = args.dataset.resolve()
    rows = json.loads((root/"catalog.json").read_text())
    plans = json.loads((root/"provenance/run_definitions.json").read_text())
    if len({(r["indent_id"], r["case"]) for r in rows}) != len(rows):
        raise ValueError("Duplicate catalog identities")
    expected = {(shape, f"p_{i:04d}") for shape, p in plans.items()
                for i in range(len(p["pressure_kPa"]))}
    if {(r["indent_id"], r["case"]) for r in rows} != expected:
        raise ValueError("Catalog does not cover the complete source pressure grid")
    def safe_path(relative):
        p = (root/relative).resolve()
        if not p.is_relative_to(root):
            raise ValueError("Catalog path leaves dataset directory")
        return p
    q_error_max, pressure_error_max, f_edge_error_max = 0., 0., 0.
    for i, row in enumerate(rows, 1):
        path = safe_path(row["gfile"])
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["gfile_sha256"]:
            raise ValueError(f"Gfile hash mismatch: {path}")
        m = json.loads(Path(str(path)+".json").read_text())
        if m["geqdsk_sha256"] != row["gfile_sha256"] or m["cocos"] != 1:
            raise ValueError("Gfile sidecar mismatch")
        d = read_geqdsk(path)
        if (check_cocos(d) != 1 or d.psirz.shape != (row["nr"], row["nz"])
                or not np.all(np.isfinite(d.psirz)) or np.any(d.pres < -1e-8)
                or np.any(d.fpol <= 0) or np.any(d.qpsi <= 0) or d.limitr != 0):
            raise ValueError(f"Invalid physical/file data: {path}")
        # Readers commonly omit a small sliver of the grid edge. Their
        # midplane inverse must still reach the actual LCFS at psi_N=1.
        psi_at = RegularGridInterpolator((d.rgrid, d.zgrid), d.psirz)
        for r in (np.linspace(d.rmaxis, d.rgrid[-1]*.999, 301),
                  np.linspace(d.rgrid[0]*1.001, d.rmaxis, 301)):
            u = (psi_at(np.column_stack((r, np.full_like(r, d.zmaxis))))-d.simag)/(d.sibry-d.simag)
            u[0 if r[0] == d.rmaxis else -1] = 0.
            mapped = interp1d(u, r, kind="cubic")(np.linspace(0, 1, 201))
            if not np.all(np.isfinite(mapped)):
                raise ValueError(f"Midplane map does not cover the LCFS: {path}")
        physics = plans[row["indent_id"]]["physics"]
        target_q = np.polynomial.polynomial.polyval(d.psin, physics["q_coefficients"])
        target_p = row["p0_kPa"]*1000*(1-d.psin**2)**2
        qe = float(np.max(abs(d.qpsi/target_q-1)))
        pe = float(np.max(abs(d.pres-target_p))/(row["p0_kPa"]*1000))
        fe = abs(float(d.fpol[-1])/m["original_metadata"]["F_edge_Tm"]-1)
        if qe > 1e-6 or pe > 1e-6 or fe > 1e-5:
            raise ValueError(f"Exported profiles disagree with source physics: {path}")
        q_error_max, pressure_error_max, f_edge_error_max = (
            max(q_error_max, qe), max(pressure_error_max, pe), max(f_edge_error_max, fe))
        if row["native_h5"]:
            native = safe_path(row["native_h5"])
            if hashlib.sha256(native.read_bytes()).hexdigest() != row["source_h5_sha256"]:
                raise ValueError("Native H5 hash mismatch")
        if geqdsk is not None:
            with gzip.open(path, "rt") as stream:
                other = geqdsk.read(stream)
            np.testing.assert_allclose(other.psi, d.psirz, rtol=0, atol=1e-14)
            for a, b in ((other.fpol, d.fpol), (other.pprime, d.pprime),
                         (other.ffprime, d.ffprim), (other.qpsi, d.qpsi)):
                np.testing.assert_allclose(a, b, rtol=0, atol=1e-14)
        if i % 64 == 0:
            print(f"Validated {i}/{len(rows)}", flush=True)
    report = dict(schema=1, cases=len(rows), families=len(plans), all_hashes_valid=True,
        all_cocos_1=True, all_midplane_maps_reach_lcfs=True,
        independent_reader="FreeQDSK" if geqdsk else None,
        max_abs_Ip_relative_difference_from_DESC=max(abs(r["Ip_relative_difference_from_DESC"]) for r in rows),
        max_F_surface_relative_std=max(r["F_surface_relative_std_max"] for r in rows),
        max_q_relative_error_from_prescribed_profile=q_error_max,
        max_pressure_error_divided_by_axis_pressure=pressure_error_max,
        max_F_edge_relative_difference_from_DESC=f_edge_error_max,
        max_inverse_mapping_residual_m=max(r["inverse_map_max_residual_m"] for r in rows),
        note="Interchange validation; does not certify a vacuum exterior or ballooning convergence")
    (root/"validation.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
