"""Standard (a48,3i4)/(5e16.9) G-EQDSK, including gzip transport.

Arrays have shape (nR, nZ); Fortran serialization makes R the fast index.
The writer preserves supplied signs and flux units. Export adapters use COCOS 1.
"""
import gzip
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np
from scipy.spatial import cKDTree

from eq_process import _parse_geqdsk, EquilibriumData


def read_geqdsk(filename):
    path = Path(filename)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        return _parse_geqdsk(stream.read(), str(path))


def _real(value, name, shape=None):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real")
    a = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(a)) or (shape is not None and a.shape != shape):
        raise ValueError(f"Invalid finite shape for {name}: expected {shape}")
    return a


def make_data(rgrid, zgrid, psi, *, axis, span, rcentr, bcentr, current,
              fpol, pres, ffprim, pprime, qpsi, boundary, header="Python GS"):
    r, z = np.asarray(rgrid), np.asarray(zgrid)
    b = np.asarray(boundary)
    return EquilibriumData(
        filename="", header=header, idum=0, nw=len(r), nh=len(z),
        rdim=float(r[-1]-r[0]), zdim=float(z[-1]-z[0]),
        rcentr=float(rcentr), rleft=float(r[0]), zmid=float((z[-1]+z[0])/2),
        rmaxis=float(axis[0]), zmaxis=float(axis[1]), simag=0., sibry=float(span),
        bcentr=float(bcentr), current=float(current),
        fpol=fpol, pres=pres, ffprim=ffprim, pprime=pprime, qpsi=qpsi,
        psirz=psi, nbbbs=len(b), limitr=0, rbbbs=b[:, 0], zbbbs=b[:, 1],
        rlim=np.empty(0), zlim=np.empty(0), rgrid=r, zgrid=z,
        psin=np.linspace(0, 1, len(r)), scalar20_raw=np.empty(0))


def write_geqdsk(data, filename, *, overwrite=False, metadata=None):
    """Write finite, uniform-grid data and, optionally, a provenance sidecar.

    Existing files are refused by default. Validation happens before writing.
    Each file is atomically replaced; a sidecar is named <gfile>.json.
    """
    path = Path(filename)
    sidecar = Path(str(path)+".json")
    if not overwrite and (path.exists() or (metadata is not None and sidecar.exists())):
        raise FileExistsError(path)
    counts = [data.nw, data.nh, data.nbbbs, data.limitr, data.idum]
    if any(isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) for n in counts):
        raise ValueError("Grid and boundary counts must be integers")
    if not (4 <= data.nw <= 999 and 4 <= data.nh <= 999
            and 3 <= data.nbbbs <= 99999 and 0 <= data.limitr <= 99999
            and -999 <= data.idum <= 9999):
        raise ValueError("Unsupported G-EQDSK grid/count bounds")
    r, z = (_real(getattr(data, k), k, (n,)) for k, n in
            (("rgrid", data.nw), ("zgrid", data.nh)))
    if r[0] <= 0 or any(np.any(np.diff(a) <= 0) or not np.allclose(
            a, np.linspace(a[0], a[-1], len(a)), rtol=0, atol=1e-12) for a in (r, z)):
        raise ValueError("G-EQDSK requires uniform, increasing, positive-R grids")
    scalar_names = ("rdim", "zdim", "rcentr", "rleft", "zmid", "rmaxis",
                    "zmaxis", "simag", "sibry", "bcentr", "current")
    scalars = {k: float(_real(getattr(data, k), k, ())) for k in scalar_names}
    if (scalars["sibry"] == scalars["simag"] or scalars["rcentr"] <= 0
            or not r[0] <= scalars["rmaxis"] <= r[-1]
            or not z[0] <= scalars["zmaxis"] <= z[-1]):
        raise ValueError("Invalid flux span or magnetic-axis/reference position")
    expected = (r[-1]-r[0], z[-1]-z[0], r[0], (z[0]+z[-1])/2)
    if not np.allclose([scalars[k] for k in ("rdim", "zdim", "rleft", "zmid")],
                       expected, rtol=0, atol=1e-12):
        raise ValueError("Grid geometry disagrees with G-EQDSK scalars")
    arrays = [_real(getattr(data, k), k, (data.nw,))
              for k in ("fpol", "pres", "ffprim", "pprime")]
    arrays += [_real(data.psirz, "psirz", (data.nw, data.nh)).ravel(order="F"),
               _real(data.qpsi, "qpsi", (data.nw,))]
    boundary = np.column_stack([_real(getattr(data, k), k, (data.nbbbs,))
                                for k in ("rbbbs", "zbbbs")])
    limiter = np.column_stack([_real(getattr(data, k), k, (data.limitr,))
                               for k in ("rlim", "zlim")])
    if np.any(boundary[:, 0] <= 0) or np.any(limiter[:, 0] <= 0):
        raise ValueError("Boundary/limiter radii must be positive")
    header = str(data.header).encode("ascii", "replace").decode().replace("\n", " ").replace("\r", " ")[:48]
    lines = [f"{header:48s}{data.idum:4d}{data.nw:4d}{data.nh:4d}\n"]
    def numbers(values):
        v = np.asarray(values).ravel()
        for first in range(0, len(v), 5):
            fields = [f"{x:16.9E}" for x in v[first:first+5]]
            if any(len(field) != 16 for field in fields):
                raise ValueError("Value exceeds G-EQDSK E16.9 format")
            lines.append("".join(fields)+"\n")
    numbers([data.rdim, data.zdim, data.rcentr, data.rleft, data.zmid,
             data.rmaxis, data.zmaxis, data.simag, data.sibry, data.bcentr,
             data.current, data.simag, 0., data.rmaxis, 0.,
             data.zmaxis, 0., data.sibry, 0., 0.])
    for a in arrays:
        numbers(a)
    lines.append(f"{data.nbbbs:5d}{data.limitr:5d}\n")
    numbers(boundary)
    numbers(limiter)
    payload = "".join(lines).encode("ascii")
    output = gzip.compress(payload, mtime=0) if path.suffix == ".gz" else payload
    extra = None
    if metadata is not None:
        extra = dict(metadata, geqdsk_sha256=hashlib.sha256(output).hexdigest())
        extra = (json.dumps(extra, indent=2, allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    def atomic_write(target, content):
        fd, name = tempfile.mkstemp(prefix=".geqdsk-", dir=target.parent)
        try:
            with os.fdopen(fd, "wb") as stream:
                stream.write(content)
            # link is an atomic no-clobber publication when overwrite=False.
            if overwrite:
                os.replace(name, target)
            elif os.name == "nt":
                # Windows rename is atomic and refuses an existing destination.
                os.rename(name, target)
            else:
                os.link(name, target)
        finally:
            if os.path.exists(name):
                os.unlink(name)
    atomic_write(path, output)
    if extra is not None:
        atomic_write(sidecar, extra)
    return path


def fill_exterior(points, inside, values, boundary, gradient, boundary_flux):
    """Normal linear continuation, for storage ONLY; this is not a vacuum solve.

    Interior samples are untouched. Positive outward flux and an explicit
    sidecar prevent a constant exterior from introducing spurious LCFS plateaus.
    """
    out = np.asarray(values).copy()
    outside = ~inside
    if not np.any(outside):
        return out
    b = np.asarray(boundary)
    grad = np.asarray(gradient)
    end = np.roll(b, -1, axis=0)
    edge = end-b
    length2 = np.sum(edge**2, axis=1)
    if np.any(length2 <= 0) or not np.all(np.isfinite(grad)):
        raise ValueError("Invalid boundary or boundary field for exterior continuation")
    queries = np.asarray(points)[outside]
    _, candidates = cKDTree((b+end)/2).query(queries, k=min(8, len(b)))
    if candidates.ndim == 1:
        candidates = candidates[:, None]
    delta = queries[:, None, :]-b[candidates]
    t = np.clip(np.sum(delta*edge[candidates], axis=2)/length2[candidates], 0, 1)
    nearest = b[candidates]+t[..., None]*edge[candidates]
    distance2 = np.sum((queries[:, None, :]-nearest)**2, axis=2)
    chosen = np.argmin(distance2, axis=1)
    rows = np.arange(len(queries))
    index, frac = candidates[rows, chosen], t[rows, chosen]
    slopes = np.linalg.norm((1-frac[:, None])*grad[index]
                           +frac[:, None]*grad[(index+1) % len(b)], axis=1)
    out[outside] = boundary_flux+slopes*np.sqrt(distance2[rows, chosen])
    if not np.all(np.isfinite(out)):
        raise ValueError("Nonfinite export grid")
    return out
