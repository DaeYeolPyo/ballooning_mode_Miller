"""GEQDSK preprocessing for local-equilibrium and s-alpha calculations."""

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import re

import contourpy
import numpy as np
from scipy.interpolate import (
    PchipInterpolator,
    RectBivariateSpline,
    interp1d,
)


# ----------------------------------------------------------------------
# Data containers
# ----------------------------------------------------------------------

@dataclass
class EquilibriumData:
    filename: str
    header: str
    idum: int
    nw: int
    nh: int

    rdim: float
    zdim: float
    rcentr: float
    rleft: float
    zmid: float

    rmaxis: float
    zmaxis: float
    simag: float
    sibry: float
    bcentr: float
    current: float

    fpol: np.ndarray
    pres: np.ndarray
    ffprim: np.ndarray
    pprime: np.ndarray
    qpsi: np.ndarray
    psirz: np.ndarray

    nbbbs: int
    limitr: int
    rbbbs: np.ndarray
    zbbbs: np.ndarray
    rlim: np.ndarray
    zlim: np.ndarray

    rgrid: np.ndarray
    zgrid: np.ndarray
    psin: np.ndarray
    scalar20_raw: np.ndarray

    @property
    def psi_axis(self):
        return self.simag

    @property
    def psi_bdry(self):
        return self.sibry


@dataclass
class FluxSurface:
    psiN: float
    psi: float

    R: np.ndarray
    Z: np.ndarray
    s_arc: np.ndarray
    dl: np.ndarray
    L: float

    F: float
    q: float
    pprime: float
    qprime: float

    Raxis: float
    Zaxis: float


@dataclass
class SurfaceIntegrals:
    psiN: float
    psi: float

    q_int: float
    q_geqdsk: float
    q_relerr: float

    V: float
    Vprime: float

    s: float
    s_hat: float
    alpha: float

    area_poloidal: float
    mean_R: float


# ----------------------------------------------------------------------
# GEQDSK reader
# ----------------------------------------------------------------------

_NUMBER_RE = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
)


def _parse_geqdsk(text, filename=""):
    lines = text.splitlines()

    if not lines:
        raise ValueError("The GEQDSK file is empty.")

    tokens = lines[0].split()

    if len(tokens) < 3:
        raise ValueError("Cannot parse the GEQDSK header.")

    try:
        idum, nw, nh = map(int, tokens[-3:])
    except ValueError as exc:
        raise ValueError(
            "Invalid idum/nw/nh in the GEQDSK header."
        ) from exc

    if nw < 2 or nh < 2:
        raise ValueError(
            "The GEQDSK grid must have at least 2 x 2 points."
        )

    raw = "\n".join(lines[1:])

    # Accept packed Fortran numbers and an optional textual diagnostic tail.
    matches = _NUMBER_RE.finditer(raw)
    position = 0
    values = []

    for match in matches:
        if raw[position:match.start()].strip():
            break

        token = match.group().replace("D", "E").replace("d", "E")
        values.append(float(token))
        position = match.end()

    if not values:
        raise ValueError("The GEQDSK file has no numeric data.")

    numbers = np.asarray(values)

    if not np.all(np.isfinite(numbers)):
        raise ValueError("GEQDSK numeric data must be finite.")

    cursor = 0

    def take(count):
        nonlocal cursor

        if cursor + count > numbers.size:
            raise ValueError(
                "Unexpected end of GEQDSK numeric data."
            )

        result = numbers[cursor:cursor + count].copy()
        cursor += count
        return result

    scalar20 = take(20)

    scalar_names = (
        "rdim", "zdim", "rcentr", "rleft", "zmid",
        "rmaxis", "zmaxis", "simag", "sibry", "bcentr",
        "current",
    )
    scalars = dict(zip(scalar_names, scalar20[:11]))

    fpol, pres, ffprim, pprime = (
        take(nw) for _ in range(4)
    )

    # MATLAB reshape(..., [nw, nh]) uses column-major ordering.
    psirz = take(nw * nh).reshape(
        (nw, nh),
        order="F",
    )

    qpsi = take(nw)
    counts = take(2)

    if (
        np.any(counts < 0)
        or np.any(counts != np.floor(counts))
    ):
        raise ValueError(
            "Boundary and limiter counts must be "
            "nonnegative integers."
        )

    nbbbs, limitr = map(int, counts)

    boundary = take(2 * nbbbs).reshape((-1, 2))
    limiter = take(2 * limitr).reshape((-1, 2))

    rgrid = np.linspace(
        scalars["rleft"],
        scalars["rleft"] + scalars["rdim"],
        nw,
    )
    zgrid = np.linspace(
        scalars["zmid"] - scalars["zdim"] / 2,
        scalars["zmid"] + scalars["zdim"] / 2,
        nh,
    )

    if scalars["sibry"] != scalars["simag"]:
        psin = np.linspace(0.0, 1.0, nw)
    else:
        psin = np.full(nw, np.nan)

    return EquilibriumData(
        filename=str(filename),
        header=" ".join(tokens[:-3]),
        idum=idum,
        nw=nw,
        nh=nh,
        **scalars,
        fpol=fpol,
        pres=pres,
        ffprim=ffprim,
        pprime=pprime,
        qpsi=qpsi,
        psirz=psirz,
        nbbbs=nbbbs,
        limitr=limitr,
        rbbbs=boundary[:, 0],
        zbbbs=boundary[:, 1],
        rlim=limiter[:, 0],
        zlim=limiter[:, 1],
        rgrid=rgrid,
        zgrid=zgrid,
        psin=psin,
        scalar20_raw=scalar20,
    )


def read_geqdsk(filename):
    """Read a GEQDSK file (plain or .gz); psirz has shape (nR, nZ)."""
    path = Path(filename)
    import gzip
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        return _parse_geqdsk(stream.read(), filename=path)


# ----------------------------------------------------------------------
# Interpolation
# ----------------------------------------------------------------------

class _Spline2D:
    """Pointwise cubic interpolation, returning NaN outside the grid."""

    def __init__(self, R, Z, values):
        self._bounds = (
            R[0], R[-1],
            Z[0], Z[-1],
        )

        self._spline = RectBivariateSpline(
            R, Z, values,
            kx=3,
            ky=3,
            s=0,
        )

    def __call__(self, R, Z):
        R, Z = np.broadcast_arrays(
            np.asarray(R, dtype=float),
            np.asarray(Z, dtype=float),
        )

        shape = R.shape
        R = R.ravel()
        Z = Z.ravel()

        rmin, rmax, zmin, zmax = self._bounds

        valid = (
            np.isfinite(R)
            & np.isfinite(Z)
            & (R >= rmin)
            & (R <= rmax)
            & (Z >= zmin)
            & (Z <= zmax)
        )

        result = np.full(R.size, np.nan)
        result[valid] = self._spline.ev(
            R[valid],
            Z[valid],
        )
        result = result.reshape(shape)

        return result.item() if result.ndim == 0 else result


class EquilibriumFields:
    """Interpolation objects with the MATLAB eqfunc interface."""

    def __init__(self, data):
        R = np.asarray(
            data.rgrid, dtype=float
        ).reshape(-1)

        Z = np.asarray(
            data.zgrid, dtype=float
        ).reshape(-1)

        psi = np.asarray(data.psirz, dtype=float)

        for grid in (R, Z):
            if (
                grid.size < 4
                or not np.all(np.isfinite(grid))
                or np.any(np.diff(grid) <= 0)
            ):
                raise ValueError(
                    "Spline grids need at least four "
                    "increasing finite points."
                )

        if (
            psi.shape != (R.size, Z.size)
            or not np.all(np.isfinite(psi))
        ):
            raise ValueError(
                "psirz must be finite and have shape (nR, nZ)."
            )

        self._psi_axis = float(data.simag)
        self._psi_span = float(data.sibry - data.simag)

        if (
            not np.isfinite(self._psi_span)
            or self._psi_span == 0
        ):
            raise ValueError(
                "sibry - simag must be finite and nonzero."
            )

        # Differentiate on the grid first, then interpolate,
        # matching the MATLAB implementation.
        psiR_grid, psiZ_grid = np.gradient(
            psi,
            R,
            Z,
            edge_order=1,
        )

        self.psi = _Spline2D(R, Z, psi)
        self.psiR = _Spline2D(R, Z, psiR_grid)
        self.psiZ = _Spline2D(R, Z, psiZ_grid)

        psiN = np.linspace(0.0, 1.0, data.nw)
        profiles = {}

        for name in (
            "fpol", "pres", "ffprim", "pprime", "qpsi"
        ):
            values = np.asarray(
                getattr(data, name),
                dtype=float,
            ).reshape(-1)

            if (
                values.size != data.nw
                or not np.all(np.isfinite(values))
            ):
                raise ValueError(
                    f"{name} must have nw finite entries."
                )

            profiles[name] = values

        # MATLAB griddedInterpolant defaults to linear interpolation.
        for output, source in (
            ("F", "fpol"),
            ("p", "pres"),
            ("FFprime", "ffprim"),
            ("pprime", "pprime"),
        ):
            interpolant = interp1d(
                psiN,
                profiles[source],
                kind="linear",
                bounds_error=False,
                fill_value="extrapolate",
                assume_sorted=True,
            )
            setattr(self, output, interpolant)

        self.q = PchipInterpolator(
            psiN,
            profiles["qpsi"],
            extrapolate=True,
        )
        self.qprimeN = self.q.derivative()

    def psiN_from_psi(self, psi):
        return (
            np.asarray(psi) - self._psi_axis
        ) / self._psi_span

    def psi_from_psiN(self, psiN):
        return (
            self._psi_axis
            + np.asarray(psiN) * self._psi_span
        )

    def qprime(self, psiN):
        """dq/dpsi, rather than dq/dpsiN."""
        return self.qprimeN(psiN) / self._psi_span

    def eval(self, R, Z):
        """Evaluate equilibrium quantities at paired or broadcast R, Z."""
        R, Z = np.broadcast_arrays(
            np.asarray(R, dtype=float),
            np.asarray(Z, dtype=float),
        )

        psi = self.psi(R, Z)
        psiR = self.psiR(R, Z)
        psiZ = self.psiZ(R, Z)

        psiN = self.psiN_from_psi(psi)
        gradPsi = np.hypot(psiR, psiZ)
        F = self.F(psiN)

        with np.errstate(
            divide="ignore",
            invalid="ignore",
        ):
            Bp = gradPsi / R
            Bphi = F / R

        return SimpleNamespace(
            R=R,
            Z=Z,
            psi=psi,
            psiR=psiR,
            psiZ=psiZ,
            psiN=psiN,
            gradPsi=gradPsi,
            F=F,
            q=self.q(psiN),
            p=self.p(psiN),
            pprime=self.pprime(psiN),
            ffprime=self.FFprime(psiN),
            qprime=self.qprime(psiN),
            Bp=Bp,
            Bphi=Bphi,
            B2=Bp**2 + Bphi**2,
        )


def build_interpolants(data):
    return EquilibriumFields(data)


# ----------------------------------------------------------------------
# Contour geometry
# ----------------------------------------------------------------------

def _polygon_area_and_rmoment(R, Z):
    Rp = np.roll(R, -1)
    Zp = np.roll(Z, -1)

    cross = R * Zp - Rp * Z

    area = 0.5 * np.sum(cross)
    r_moment = np.sum((R + Rp) * cross) / 6.0

    return area, r_moment


def _remove_duplicate_endpoint(points):
    scale = max(
        1.0,
        np.max(np.hypot(points[:, 0], points[:, 1])),
    )

    if (
        points.shape[0] > 1
        and np.linalg.norm(points[0] - points[-1])
        < 1e-10 * scale
    ):
        points = points[:-1]

    return points


def _contains_point(points, x, y):
    """Ray-crossing polygon test, including the boundary."""
    a = points
    b = np.roll(points, -1, axis=0)

    edge = b - a
    relative = np.array([x, y]) - a

    tolerance = 1e-12 * max(
        1.0,
        np.max(np.abs(points)),
    )

    cross = (
        edge[:, 0] * relative[:, 1]
        - edge[:, 1] * relative[:, 0]
    )

    on_edge = (
        (
            np.abs(cross)
            <= tolerance * np.linalg.norm(edge, axis=1)
        )
        & (
            np.sum(
                relative * (np.array([x, y]) - b),
                axis=1,
            )
            <= tolerance**2
        )
    )

    if np.any(on_edge):
        return True

    crossing = (
        (a[:, 1] > y) != (b[:, 1] > y)
    )

    xcross = (
        a[crossing, 0]
        + (y - a[crossing, 1])
        * edge[crossing, 0]
        / edge[crossing, 1]
    )

    return bool(np.count_nonzero(xcross > x) % 2)


def _orient_surface(points, orientation):
    if orientation not in (
        "ccw", "cw", "outboard_up"
    ):
        raise ValueError(
            "orientation must be 'ccw', 'cw', "
            "or 'outboard_up'."
        )

    area, _ = _polygon_area_and_rmoment(
        points[:, 0],
        points[:, 1],
    )

    if (area > 0) != (orientation != "cw"):
        points = points[::-1]

    return points


def _start_surface_at_outboard(
    points,
    Zaxis,
    force_up=False,
):
    R, Z = points.T
    eps = np.finfo(float).eps

    score = (
        ((R - R.max()) / max(np.ptp(R), eps))**2
        + ((Z - Zaxis) / max(np.ptp(Z), eps))**2
    )

    points = np.roll(
        points,
        -int(np.argmin(score)),
        axis=0,
    )

    if force_up and points[1, 1] < points[0, 1]:
        points = np.concatenate(
            (points[:1], points[:0:-1])
        )

    return points


def _resample_closed_curve(points, npoints):
    closed = np.vstack((points, points[0]))

    s = np.r_[
        0.0,
        np.cumsum(
            np.linalg.norm(
                np.diff(closed, axis=0),
                axis=1,
            )
        ),
    ]

    if not np.isfinite(s[-1]) or s[-1] <= 0:
        raise ValueError(
            "The contour has zero or invalid length."
        )

    keep = np.r_[True, np.diff(s) > 0]
    s = s[keep]
    closed = closed[keep]

    sq = np.linspace(
        0.0,
        s[-1],
        npoints,
        endpoint=False,
    )

    R = PchipInterpolator(
        s, closed[:, 0]
    )(sq)

    Z = PchipInterpolator(
        s, closed[:, 1]
    )(sq)

    # Includes the final segment from the last point to the first.
    dl = np.hypot(
        np.roll(R, -1) - R,
        np.roll(Z, -1) - Z,
    )

    return R, Z, sq, dl, float(s[-1])


# ----------------------------------------------------------------------
# Flux surfaces
# ----------------------------------------------------------------------

def extract_flux_surface(
    data,
    fields,
    target_psiN,
    *,
    npoints=1024,
    start="outboard",
    orientation="outboard_up",
):
    if (
        not isinstance(npoints, (int, np.integer))
        or npoints < 4
    ):
        raise ValueError(
            "npoints must be an integer >= 4."
        )

    if start not in ("outboard", "native"):
        raise ValueError(
            "start must be 'outboard' or 'native'."
        )

    if orientation not in (
        "ccw", "cw", "outboard_up"
    ):
        raise ValueError("Unknown surface orientation.")

    target_psiN = float(target_psiN)

    if not np.isfinite(target_psiN):
        raise ValueError("target_psiN must be finite.")

    target_psi = (
        data.simag
        + target_psiN * (data.sibry - data.simag)
    )

    # ContourPy expects z.shape == (nZ, nR).
    generator = contourpy.contour_generator(
        x=data.rgrid,
        y=data.zgrid,
        z=np.ma.masked_invalid(data.psirz.T),
        name="serial",
        line_type="Separate",
        corner_mask=False,
    )

    best = None
    best_area = -np.inf

    for points in generator.lines(target_psi):
        if (
            len(points) < 4
            or not np.all(np.isfinite(points))
        ):
            continue

        scale = max(
            1.0,
            np.max(
                np.hypot(points[:, 0], points[:, 1])
            ),
        )

        # Test closure before removing the repeated endpoint.
        if (
            np.linalg.norm(points[0] - points[-1])
            >= 1e-8 * scale
        ):
            continue

        points = _remove_duplicate_endpoint(points)

        if (
            len(points) < 3
            or not _contains_point(
                points,
                data.rmaxis,
                data.zmaxis,
            )
        ):
            continue

        area, _ = _polygon_area_and_rmoment(
            points[:, 0],
            points[:, 1],
        )

        if abs(area) > best_area:
            best = points
            best_area = abs(area)

    if best is None:
        raise ValueError(
            "No closed contour enclosing the magnetic "
            f"axis at psiN={target_psiN:g}."
        )

    scale = max(
        1.0,
        np.max(np.hypot(best[:, 0], best[:, 1])),
    )

    keep = np.r_[
        True,
        np.linalg.norm(
            np.diff(best, axis=0),
            axis=1,
        ) > 1e-10 * scale,
    ]

    best = _remove_duplicate_endpoint(best[keep])

    if len(best) < 3:
        raise ValueError(
            "Too few distinct contour points."
        )

    best = _orient_surface(best, orientation)

    if start == "outboard":
        best = _start_surface_at_outboard(
            best,
            data.zmaxis,
            force_up=(orientation == "outboard_up"),
        )

    R, Z, s_arc, dl, L = _resample_closed_curve(
        best,
        npoints,
    )

    return FluxSurface(
        psiN=target_psiN,
        psi=float(target_psi),
        R=R,
        Z=Z,
        s_arc=s_arc,
        dl=dl,
        L=L,
        F=float(fields.F(target_psiN)),
        q=float(fields.q(target_psiN)),
        pprime=float(fields.pprime(target_psiN)),
        qprime=float(fields.qprime(target_psiN)),
        Raxis=float(data.rmaxis),
        Zaxis=float(data.zmaxis),
    )


# ----------------------------------------------------------------------
# Contour integrals
# ----------------------------------------------------------------------

def compute_contour_integrals(
    data,
    fields,
    surface,
    *,
    is_normalized=False,
):
    """
    Keep the MATLAB quadrature and positive Vprime convention.

    The data argument is retained for compatibility with the MATLAB API.
    """
    R = surface.R
    Z = surface.Z
    dl = surface.dl

    values = fields.eval(R, Z)
    gradPsi = values.gradPsi

    if (
        np.any(~np.isfinite(gradPsi))
        or np.any(gradPsi <= 0)
    ):
        raise ValueError(
            "Invalid |grad psi| on the flux surface."
        )

    if (
        np.any(R <= 0)
        or np.any(~np.isfinite(values.F))
    ):
        raise ValueError(
            "R must be positive and F must be finite."
        )

    if (
        surface.q == 0
        or not np.isfinite(surface.q)
        or surface.Raxis <= 0
    ):
        raise ValueError(
            "Surface q and Raxis are invalid."
        )

    q_int = float(
        np.sum(
            values.F / (R * gradPsi) * dl
        ) / (2 * np.pi)
    )

    q_geqdsk = float(fields.q(surface.psiN))

    Vprime = float(
        2 * np.pi * np.sum(R / gradPsi * dl)
    )

    area, r_moment = _polygon_area_and_rmoment(R, Z)
    V = float(2 * np.pi * abs(r_moment))

    if (
        not np.isfinite(V + Vprime)
        or V <= 0
        or Vprime <= 0
    ):
        raise ValueError(
            "V and Vprime must be positive and finite."
        )

    mu0 = (
        1.0
        if is_normalized
        else 4 * np.pi * 1e-7
    )

    s = 2 * V * surface.qprime / Vprime

    alpha = (
        -(2 * Vprime / (4 * np.pi**2))
        * np.sqrt(
            V / (2 * np.pi**2 * surface.Raxis)
        )
        * mu0
        * surface.pprime
    )

    return SurfaceIntegrals(
        psiN=surface.psiN,
        psi=surface.psi,
        q_int=q_int,
        q_geqdsk=q_geqdsk,
        q_relerr=(
            (q_int - q_geqdsk)
            / max(1.0, abs(q_geqdsk))
        ),
        V=V,
        Vprime=Vprime,
        s=float(s),
        s_hat=float(s / surface.q),
        alpha=float(alpha),
        area_poloidal=float(abs(area)),
        mean_R=float(
            abs(r_moment)
            / max(abs(area), np.finfo(float).eps)
        ),
    )


def s_to_qprime(s, integrals):
    """Here s = q*s_hat, matching the original MATLAB function."""
    return (
        np.asarray(s)
        * integrals.Vprime
        / (2 * integrals.V)
    )


# ----------------------------------------------------------------------
# COCOS sign check
# ----------------------------------------------------------------------

def check_cocos(data, sigma_RphiZ=1, e_Bp=0):
    """
    Port the original sign check.

    This function does not convert the equilibrium convention.
    """
    if (
        sigma_RphiZ not in (-1, 1)
        or e_Bp not in (0, 1)
    ):
        raise ValueError(
            "sigma_RphiZ must be +/-1 "
            "and e_Bp must be 0 or 1."
        )

    def main_sign(values):
        values = np.asarray(values)
        values = values[np.isfinite(values)]

        if values.size == 0:
            return np.nan

        tolerance = 1e-10 * max(
            1.0,
            np.max(np.abs(values)),
        )

        values = values[
            np.abs(values) > tolerance
        ]

        if values.size == 0:
            return np.nan

        return np.sign(np.median(values))

    sigma_Ip = np.sign(data.current)
    sigma_B0 = np.sign(data.bcentr)

    sigma_Bp = (
        sigma_Ip
        * np.sign(data.sibry - data.simag)
    )

    sigma_rhotp = (
        sigma_Ip
        * sigma_B0
        * main_sign(data.qpsi)
    )

    if (
        main_sign(data.fpol) != sigma_B0
        or main_sign(data.pprime)
        != -sigma_Ip * sigma_Bp
    ):
        raise ValueError(
            "COCOS convention signs are not valid."
        )

    base = {
        (1, 1): 1,
        (-1, -1): 3,
        (1, -1): 5,
        (-1, 1): 7,
    }.get(
        (sigma_Bp, sigma_rhotp),
        np.nan,
    )

    return (
        base
        + int(sigma_RphiZ == -1)
        + 10 * e_Bp
    )


# ----------------------------------------------------------------------
# User-facing interface
# ----------------------------------------------------------------------

class Equilibrium:
    """Reusable equilibrium data and field evaluator."""

    def __init__(self, data, fields=None):
        self.data = data

        # A future FEM adapter can supply its own field evaluator.
        self.fields = (
            build_interpolants(data)
            if fields is None
            else fields
        )

    @classmethod
    def from_geqdsk(cls, filename):
        return cls(read_geqdsk(filename))

    def extract_flux_surface(
        self,
        psiN,
        *,
        npoints=1024,
        start="outboard",
        orientation="outboard_up",
    ):
        return extract_flux_surface(
            self.data,
            self.fields,
            psiN,
            npoints=npoints,
            start=start,
            orientation=orientation,
        )

    def compute_contour_integrals(
        self,
        surface,
        *,
        is_normalized=False,
    ):
        return compute_contour_integrals(
            self.data,
            self.fields,
            surface,
            is_normalized=is_normalized,
        )


# ----------------------------------------------------------------------
# Standalone execution
# ----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument("geqdsk", type=Path)
    parser.add_argument(
        "--psiN",
        type=float,
        default=0.7,
    )
    parser.add_argument(
        "--npoints",
        type=int,
        default=1024,
    )

    args = parser.parse_args()

    eq = Equilibrium.from_geqdsk(args.geqdsk)

    surface = eq.extract_flux_surface(
        args.psiN,
        npoints=args.npoints,
    )
    integrals = eq.compute_contour_integrals(surface)

    print(f"psiN = {surface.psiN:.6g}")
    print(
        f"q(profile) = {surface.q:.9g}, "
        f"q(integral) = {integrals.q_int:.9g}"
    )
    print(
        f"V = {integrals.V:.9g}, "
        f"Vprime = {integrals.Vprime:.9g}"
    )
    print(
        f"s_hat = {integrals.s_hat:.9g}, "
        f"alpha = {integrals.alpha:.9g}"
    )
