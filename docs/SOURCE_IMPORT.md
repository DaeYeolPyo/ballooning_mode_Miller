# Source provenance

Imported on 2026-11-19 from the user-supplied server source archive:

```text
python_sources_for_github.tar.gz
SHA-256: 96a1098cc9390d180be465f228040bc8bcf4c22c9d21b575c55f8b2a4f65916c
```

Retained: `gs_fixed`, `eq_process.py`, `straight_field_line.py`,
`local_equilibrium.py`, `salpha.py`, `sturm_liouville.py`, their core tests and
two GS examples. The computational FEM, mesh, profiles, Picard/prescribed-q,
local-equilibrium and P2 eigenvalue formulations are unchanged.

New: gfile writer/export adapters, archive conversion/catalog CLI, handoff
dataset, examples, packaging and validation. `GSEquilibrium.export_gfile` is a
thin method; `eq_process.read_geqdsk` additionally accepts gzip. Native DESC
H5 decoding comes from `examples/chance_h5_geometry.py` in the same archive;
its coefficient/mode/energy functions are retained in `ballooning/desc.py`.
Experimental coarse-mapping scans, scheduler jobs, result dumps and MATLAB
source are omitted from the current tree. The MATLAB Git history is retained.

For the declared NumPy >=1.24 support, `local_equilibrium.py` and `salpha.py`
use SciPy's equivalent `trapezoid` function in place of the NumPy-2-only
`np.trapezoid`. The quadrature formula is unchanged.

The bean dataset originates from the complete Chance reproduction run with
L64/M40/N0 DESC equilibria. Per-file H5 hashes are checked against their
original metadata. Physical input definitions and the original scan summaries
are preserved under `data/bean/provenance`. No new equilibrium or whole-plasma
stability scan was run to create the gfiles.
