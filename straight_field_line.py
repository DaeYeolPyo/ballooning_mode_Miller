"""PEST, Hamada, and Boozer coordinates and metric calculations."""

from dataclasses import dataclass
from types import SimpleNamespace
import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator

from eq_process import FluxSurface, SurfaceIntegrals


# ----------------------------------------------------------------------
# Data containers
# ----------------------------------------------------------------------

@dataclass
class SFLAngles:
    psi: float
    psiN: float

    theta_geom: np.ndarray
    theta_PEST: np.ndarray
    theta_Hamada: np.ndarray
    theta_Boozer: np.ndarray

    lambda_Hamada: np.ndarray
    lambda_Boozer: np.ndarray

    q_used: float
    weights: SimpleNamespace
    diagnostics: SimpleNamespace


@dataclass
class SFLSurface:
    surf: FluxSurface
    ints: SurfaceIntegrals
    ang: SFLAngles


@dataclass
class SFLMap:
    theta: np.ndarray

    R: np.ndarray
    Z: np.ndarray
    lambda_shift: np.ndarray

    psi: np.ndarray
    psiN: np.ndarray
    q: np.ndarray
    F: np.ndarray
    pprime: np.ndarray
    V: np.ndarray
    Vprime: np.ndarray


@dataclass
class SFLMaps:
    PEST: SFLMap
    Hamada: SFLMap
    Boozer: SFLMap


@dataclass
class SymmetricTensor:
    psipsi: np.ndarray
    psitheta: np.ndarray
    psizeta: np.ndarray
    thetatheta: np.ndarray
    thetazeta: np.ndarray
    zetazeta: np.ndarray


@dataclass
class VectorComponents:
    psi: np.ndarray
    theta: np.ndarray
    zeta: np.ndarray


@dataclass
class MetricData:
    Jacobian: np.ndarray
    Jacobian_abs: np.ndarray
    detg: np.ndarray

    g_cov: SymmetricTensor
    g_contra: SymmetricTensor

    B_contra: VectorComponents
    B_cov: VectorComponents
    B2: np.ndarray
    B2_direct: np.ndarray
    B2_metric: np.ndarray

    J_contra: VectorComponents
    J_cov: VectorComponents


# ----------------------------------------------------------------------
# Common helpers
# ----------------------------------------------------------------------

def _vector(value, name, minimum=1):
    value = np.asarray(value, dtype=float)

    if value.ndim not in (1, 2) or (
        value.ndim == 2 and 1 not in value.shape
    ):
        raise ValueError(f"{name} must be a vector.")

    value = value.reshape(-1)

    if (
        value.size < minimum
        or not np.all(np.isfinite(value))
    ):
        raise ValueError(
            f"{name} must have at least {minimum} finite entries."
        )

    return value


def _normalized_cumulative_angle(weight, dl):
    # dl[i] connects point i to point i+1, including periodic closure.
    segment_integrals = (
        0.5 * (weight + np.roll(weight, -1)) * dl
    )
    total = np.sum(segment_integrals)

    if (
        not np.isfinite(total)
        or abs(total) < np.finfo(float).eps
    ):
        raise ValueError(
            "The angle weight integral is zero or invalid."
        )

    cumulative = np.r_[
        0.0,
        np.cumsum(segment_integrals[:-1]),
    ]

    theta = np.mod(
        2 * np.pi * cumulative / total,
        2 * np.pi,
    )
    theta[0] = 0.0

    theta = np.unwrap(theta)

    return theta - theta[0]


def _geometric_poloidal_angle(R, Z, Raxis, Zaxis):
    theta = np.unwrap(
        np.arctan2(Z - Zaxis, R - Raxis)
    )
    theta = theta - theta[0]

    if theta[-1] < 0:
        theta = theta + 2 * np.pi

    theta = np.mod(theta, 2 * np.pi)
    theta[0] = 0.0

    return np.unwrap(theta)


def _unwrap_periodic_difference(value):
    value = np.unwrap(value)
    return value - value[0]


# ----------------------------------------------------------------------
# Angles on one flux surface
# ----------------------------------------------------------------------

def straight_field_line_angles(
    fields,
    surface,
    integrals,
    *,
    use_q="integral",
):
    """Compute angles on one surface; endpoints are not duplicated."""
    R = _vector(surface.R, "R", 4)
    Z = _vector(surface.Z, "Z", 4)
    dl = _vector(surface.dl, "dl", 4)

    if not (R.size == Z.size == dl.size):
        raise ValueError(
            "R, Z, and dl must have equal lengths."
        )

    if np.any(R <= 0) or np.any(dl <= 0):
        raise ValueError("R and dl must be positive.")

    values = fields.eval(R, Z)

    gradPsi, F, B2 = (
        np.broadcast_to(
            np.asarray(value, dtype=float),
            R.shape,
        )
        for value in (
            values.gradPsi,
            values.F,
            values.B2,
        )
    )

    if not all(
        np.all(np.isfinite(value))
        for value in (gradPsi, F, B2)
    ):
        raise ValueError(
            "The sampled equilibrium fields must be finite."
        )

    if np.any(gradPsi <= 0) or np.any(B2 <= 0):
        raise ValueError(
            "|grad psi| and B2 must be positive."
        )

    wPEST = F / (R * gradPsi)
    wHamada = R / gradPsi
    wBoozer = R * B2 / gradPsi

    theta_PEST = _normalized_cumulative_angle(
        wPEST, dl
    )
    theta_Hamada = _normalized_cumulative_angle(
        wHamada, dl
    )
    theta_Boozer = _normalized_cumulative_angle(
        wBoozer, dl
    )

    use_q = use_q.lower()

    if use_q == "integral":
        q = float(integrals.q_int)
    elif use_q == "geqdsk":
        q = float(integrals.q_geqdsk)
    else:
        raise ValueError(
            "use_q must be 'integral' or 'geqdsk'."
        )

    if not np.isfinite(q):
        raise ValueError(
            "The selected q must be finite."
        )

    lambda_Hamada = q * _unwrap_periodic_difference(
        theta_Hamada - theta_PEST
    )
    lambda_Boozer = q * _unwrap_periodic_difference(
        theta_Boozer - theta_PEST
    )

    return SFLAngles(
        psi=surface.psi,
        psiN=surface.psiN,
        theta_geom=_geometric_poloidal_angle(
            R,
            Z,
            surface.Raxis,
            surface.Zaxis,
        ),
        theta_PEST=theta_PEST,
        theta_Hamada=theta_Hamada,
        theta_Boozer=theta_Boozer,
        lambda_Hamada=lambda_Hamada,
        lambda_Boozer=lambda_Boozer,
        q_used=q,
        weights=SimpleNamespace(
            PEST=wPEST,
            Hamada=wHamada,
            Boozer=wBoozer,
        ),
        diagnostics=SimpleNamespace(
            PEST_monotone=bool(
                np.all(np.diff(theta_PEST) > 0)
            ),
            Hamada_monotone=bool(
                np.all(np.diff(theta_Hamada) > 0)
            ),
            Boozer_monotone=bool(
                np.all(np.diff(theta_Boozer) > 0)
            ),
        ),
    )


def get_raw_surface_data(surface, name):
    """Return theta, R, Z, lambda_shift from one SFLSurface."""
    surf = surface.surf
    ang = surface.ang
    name = name.lower()

    if name == "pest":
        theta = ang.theta_PEST
        shift = np.zeros_like(surf.R)

    elif name == "hamada":
        theta = ang.theta_Hamada
        shift = ang.lambda_Hamada

    elif name == "boozer":
        theta = ang.theta_Boozer
        shift = ang.lambda_Boozer

    else:
        raise ValueError(
            "name must be 'PEST', 'Hamada', or 'Boozer'."
        )

    return theta, surf.R, surf.Z, shift


# ----------------------------------------------------------------------
# Remapping to uniform theta
# ----------------------------------------------------------------------

def _remap_single_surface(
    R,
    Z,
    shift,
    theta_raw,
    theta_uniform,
):
    arrays = [
        _vector(value, name, 4)
        for value, name in (
            (R, "R"),
            (Z, "Z"),
            (shift, "lambda_shift"),
            (theta_raw, "theta_raw"),
        )
    ]

    R, Z, shift, theta = arrays

    if any(value.size != R.size for value in arrays):
        raise ValueError(
            "All raw surface arrays must have equal lengths."
        )

    theta = np.unwrap(theta - theta[0])
    order = np.argsort(theta, kind="stable")

    theta = theta[order]
    values = np.column_stack((R, Z, shift))[order]

    keep = np.r_[
        True,
        np.diff(theta) > 1e-12,
    ]
    theta = theta[keep]
    values = values[keep]

    if (
        theta.size < 3
        or abs(theta[0]) > 1e-12
        or theta[-1] >= 2 * np.pi
    ):
        raise ValueError(
            "Raw angles must span one turn "
            "without a duplicate endpoint."
        )

    # Remove the arbitrary constant toroidal shift.
    values[:, 2] -= values[0, 2]

    # Add the periodic endpoint only for interpolation.
    theta_ext = np.r_[theta, 2 * np.pi]
    values_ext = np.vstack((values, values[0]))
    values_ext[-1, 2] = 0.0

    mapped = PchipInterpolator(
        theta_ext,
        values_ext,
        axis=0,
    )(theta_uniform)

    return (
        mapped[:, 0],
        mapped[:, 1],
        mapped[:, 2],
    )


def build_sfl_maps(surfaces, ntheta=256):
    """Map surfaces to uniform theta.

    R, Z, and lambda_shift have shape (n_surface, ntheta).
    The order of the input surfaces is preserved.
    """
    surfaces = list(surfaces)

    if not surfaces:
        raise ValueError(
            "At least one SFLSurface is required."
        )

    if (
        not isinstance(ntheta, (int, np.integer))
        or ntheta < 4
    ):
        raise ValueError(
            "ntheta must be an integer >= 4."
        )

    theta = np.linspace(
        0.0,
        2 * np.pi,
        ntheta,
        endpoint=False,
    )

    metadata = {
        "psi": np.array(
            [item.surf.psi for item in surfaces]
        ),
        "psiN": np.array(
            [item.surf.psiN for item in surfaces]
        ),
        "q": np.array(
            [item.ang.q_used for item in surfaces]
        ),
        "F": np.array(
            [item.surf.F for item in surfaces]
        ),
        "pprime": np.array(
            [item.surf.pprime for item in surfaces]
        ),
        "V": np.array(
            [item.ints.V for item in surfaces]
        ),
        "Vprime": np.array(
            [item.ints.Vprime for item in surfaces]
        ),
    }

    if not all(
        np.all(np.isfinite(value))
        for value in metadata.values()
    ):
        raise ValueError(
            "Surface metadata must be finite."
        )

    results = {}

    for name in ("PEST", "Hamada", "Boozer"):
        Rmap = np.empty((len(surfaces), ntheta))
        Zmap = np.empty_like(Rmap)
        shifts = np.empty_like(Rmap)

        for k, surface in enumerate(surfaces):
            raw_theta, R, Z, shift = get_raw_surface_data(
                surface,
                name,
            )

            if np.any(np.diff(raw_theta) <= 0):
                raise ValueError(
                    f"{name} angles must be strictly increasing."
                )

            Rmap[k], Zmap[k], shifts[k] = (
                _remap_single_surface(
                    R,
                    Z,
                    shift,
                    raw_theta,
                    theta,
                )
            )

        results[name] = SFLMap(
            theta=theta.copy(),
            R=Rmap,
            Z=Zmap,
            lambda_shift=shifts,
            **{
                key: value.copy()
                for key, value in metadata.items()
            },
        )

    return SFLMaps(**results)


# ----------------------------------------------------------------------
# Finite differences
# ----------------------------------------------------------------------

def _d_dpsi(A, psi):
    """Match the original secant formula, including nonuniform psi."""
    result = np.empty_like(A)

    result[1:-1] = (
        (A[2:] - A[:-2])
        / (psi[2:] - psi[:-2])[:, None]
    )

    result[0] = (
        (A[1] - A[0])
        / (psi[1] - psi[0])
    )

    result[-1] = (
        (A[-1] - A[-2])
        / (psi[-1] - psi[-2])
    )

    return result


def _d_dtheta(A, theta):
    step = theta[1] - theta[0]

    return (
        np.roll(A, -1, axis=1)
        - np.roll(A, 1, axis=1)
    ) / (2 * step)


# ----------------------------------------------------------------------
# Metric, magnetic field, and current
# ----------------------------------------------------------------------

def compute_metrics(
    sfl_map,
    fields,
    *,
    mu0=4 * np.pi * 1e-7,
):
    """Compute metric, B, and J using the MATLAB conventions."""
    psi = _vector(sfl_map.psi, "psi", 2)
    theta = _vector(sfl_map.theta, "theta", 4)

    delta_psi = np.diff(psi)

    if not (
        np.all(delta_psi > 0)
        or np.all(delta_psi < 0)
    ):
        raise ValueError(
            "psi must be strictly increasing or decreasing."
        )

    delta_theta = np.diff(
        np.r_[theta, theta[0] + 2 * np.pi]
    )

    if not np.allclose(
        delta_theta,
        2 * np.pi / theta.size,
        rtol=1e-10,
        atol=1e-12,
    ):
        raise ValueError(
            "theta must be uniform over one period "
            "without a duplicate endpoint."
        )

    if not np.isfinite(mu0) or mu0 <= 0:
        raise ValueError(
            "mu0 must be positive and finite."
        )

    R, Z, lam = (
        np.asarray(value, dtype=float)
        for value in (
            sfl_map.R,
            sfl_map.Z,
            sfl_map.lambda_shift,
        )
    )

    shape = (psi.size, theta.size)

    if any(
        value.shape != shape
        for value in (R, Z, lam)
    ):
        raise ValueError(
            "R, Z, and lambda_shift must have "
            "shape (n_surface, ntheta)."
        )

    if not all(
        np.all(np.isfinite(value))
        for value in (R, Z, lam)
    ):
        raise ValueError(
            "Mapped R, Z, and lambda_shift must be finite."
        )

    if np.any(R <= 0):
        raise ValueError(
            "Mapped R must be positive."
        )

    F = _vector(sfl_map.F, "F")
    q = _vector(sfl_map.q, "q")

    if F.size != psi.size or q.size != psi.size:
        raise ValueError(
            "F and q need one value per flux surface."
        )

    F2d = F[:, None]
    q2d = q[:, None]

    psiR = np.broadcast_to(
        fields.psiR(R, Z),
        shape,
    )
    psiZ = np.broadcast_to(
        fields.psiZ(R, Z),
        shape,
    )

    if (
        not np.all(np.isfinite(psiR))
        or not np.all(np.isfinite(psiZ))
    ):
        raise ValueError(
            "The interpolated flux gradients must be finite."
        )

    R_psi, Z_psi, l_psi = (
        _d_dpsi(value, psi)
        for value in (R, Z, lam)
    )
    R_th, Z_th, l_th = (
        _d_dtheta(value, theta)
        for value in (R, Z, lam)
    )

    Jacobian = R * (
        R_psi * Z_th - Z_psi * R_th
    )

    if np.median(Jacobian) < 0:
        warnings.warn(
            "Signed Jacobian is mostly negative; "
            "its absolute value is used as in MATLAB.",
            RuntimeWarning,
            stacklevel=2,
        )

    Jacobian_abs = np.abs(Jacobian)

    if (
        np.any(Jacobian_abs == 0)
        or not np.all(np.isfinite(Jacobian_abs))
    ):
        raise ValueError(
            "The coordinate Jacobian is zero or invalid."
        )

    # Covariant metric.
    g11 = (
        R_psi**2
        + Z_psi**2
        + R**2 * l_psi**2
    )
    g12 = (
        R_psi * R_th
        + Z_psi * Z_th
        + R**2 * l_psi * l_th
    )
    g13 = -R**2 * l_psi

    g22 = (
        R_th**2
        + Z_th**2
        + R**2 * l_th**2
    )
    g23 = -R**2 * l_th
    g33 = R**2

    detg = (
        g11 * (g22 * g33 - g23**2)
        - g12 * (g12 * g33 - g13 * g23)
        + g13 * (g12 * g23 - g13 * g22)
    )

    if (
        not np.all(np.isfinite(detg))
        or np.any(
            np.abs(detg) < np.finfo(float).eps
        )
    ):
        raise ValueError(
            "The metric tensor is singular or invalid."
        )

    # Contravariant metric.
    gcon11 = (g22 * g33 - g23**2) / detg
    gcon12 = (g13 * g23 - g12 * g33) / detg
    gcon13 = (g12 * g23 - g13 * g22) / detg

    gcon22 = (g11 * g33 - g13**2) / detg
    gcon23 = (g12 * g13 - g11 * g23) / detg
    gcon33 = (g11 * g22 - g12**2) / detg

    # Magnetic field components.
    Bcon = VectorComponents(
        psi=np.zeros_like(R),
        theta=1.0 / Jacobian_abs,
        zeta=np.broadcast_to(
            q2d / Jacobian_abs,
            shape,
        ).copy(),
    )

    Bcov = VectorComponents(
        psi=(
            g12 * Bcon.theta
            + g13 * Bcon.zeta
        ),
        theta=(
            g22 * Bcon.theta
            + g23 * Bcon.zeta
        ),
        zeta=(
            g23 * Bcon.theta
            + g33 * Bcon.zeta
        ),
    )

    B2_direct = (
        psiR**2 + psiZ**2 + F2d**2
    ) / R**2

    B2_metric = (
        Bcov.psi * Bcon.psi
        + Bcov.theta * Bcon.theta
        + Bcov.zeta * Bcon.zeta
    )

    # Current density from curl B.
    dBpsi_dth = _d_dtheta(Bcov.psi, theta)
    dBtheta_dpsi = _d_dpsi(Bcov.theta, psi)
    dBzeta_dpsi = _d_dpsi(Bcov.zeta, psi)
    dBzeta_dth = _d_dtheta(Bcov.zeta, theta)

    Jcon = VectorComponents(
        psi=(
            dBzeta_dth / Jacobian_abs / mu0
        ),
        theta=(
            -dBzeta_dpsi / Jacobian_abs / mu0
        ),
        zeta=(
            (dBtheta_dpsi - dBpsi_dth)
            / Jacobian_abs
            / mu0
        ),
    )

    Jcov = VectorComponents(
        psi=(
            g11 * Jcon.psi
            + g12 * Jcon.theta
            + g13 * Jcon.zeta
        ),
        theta=(
            g12 * Jcon.psi
            + g22 * Jcon.theta
            + g23 * Jcon.zeta
        ),
        zeta=(
            g13 * Jcon.psi
            + g23 * Jcon.theta
            + g33 * Jcon.zeta
        ),
    )

    return MetricData(
        Jacobian=Jacobian,
        Jacobian_abs=Jacobian_abs,
        detg=detg,
        g_cov=SymmetricTensor(
            g11, g12, g13,
            g22, g23, g33,
        ),
        g_contra=SymmetricTensor(
            gcon11, gcon12, gcon13,
            gcon22, gcon23, gcon33,
        ),
        B_contra=Bcon,
        B_cov=Bcov,
        B2=B2_direct,
        B2_direct=B2_direct,
        B2_metric=B2_metric,
        J_contra=Jcon,
        J_cov=Jcov,
    )


# ----------------------------------------------------------------------
# Standalone execution
# ----------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    from eq_process import Equilibrium

    parser = argparse.ArgumentParser(
        description=__doc__
    )
    parser.add_argument("geqdsk")
    parser.add_argument(
        "--nsurfaces",
        type=int,
        default=25,
    )
    parser.add_argument(
        "--npoints",
        type=int,
        default=1024,
    )
    parser.add_argument(
        "--ntheta",
        type=int,
        default=256,
    )

    args = parser.parse_args()

    if args.nsurfaces < 2:
        parser.error(
            "--nsurfaces must be at least 2"
        )

    eq = Equilibrium.from_geqdsk(args.geqdsk)
    surfaces = []

    for psiN in np.linspace(
        0.1, 0.9, args.nsurfaces
    ):
        surf = eq.extract_flux_surface(
            psiN,
            npoints=args.npoints,
        )
        ints = eq.compute_contour_integrals(surf)

        ang = straight_field_line_angles(
            eq.fields,
            surf,
            ints,
        )

        surfaces.append(
            SFLSurface(surf, ints, ang)
        )

    maps = build_sfl_maps(
        surfaces,
        ntheta=args.ntheta,
    )

    for name in ("PEST", "Hamada", "Boozer"):
        metric = compute_metrics(
            getattr(maps, name),
            eq.fields,
        )

        difference = (
            metric.B2_metric - metric.B2_direct
        )
        error = (
            np.linalg.norm(difference)
            / np.linalg.norm(metric.B2_direct)
        )

        print(
            f"{name}: B2 relative RMS difference "
            f"= {error:.6e}"
        )
