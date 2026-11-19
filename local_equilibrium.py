"""Local-equilibrium fits, Mercier-Luc geometry, and ballooning coefficients.

Dependency: eq_process.py from the previous port.
Angles in Boundary/LocalEquilibrium include both 0 and 2*pi.
Functions fit_* accept (eq, psiN) or (data, fields, psiN).
Fit option dpsi is a half-width in normalized flux psiN.
Nonlinear fits expose fit.converged; max_iterations need not imply convergence.
The MATLAB scalar-curvature convention is retained. Inspect vector_consistent
before using coefficients from a geometry with a different normal orientation.
"""

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Optional
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import minimize_scalar

from eq_process import extract_flux_surface

_EPS = np.finfo(float).eps
_TWO_PI = 2 * np.pi
_MU0 = 4 * np.pi * 1e-7


@dataclass
class Boundary:
    theta: np.ndarray
    R: np.ndarray
    Z: np.ndarray
    Rt: np.ndarray
    Zt: np.ndarray
    Rtt: np.ndarray
    Ztt: np.ndarray
    Req: Optional[np.ndarray] = None
    Zeq: Optional[np.ndarray] = None


class ShapeParameters(SimpleNamespace):
    """Model-specific coefficients, radial derivatives, and fit diagnostics."""


class LocalEquilibrium(SimpleNamespace):
    """Common closed-grid geometry consumed by evaluate_theta_psi."""


class ThetaDiagnostic(SimpleNamespace):
    """Geometry and radial equilibrium on the closed PEST grid."""


class BallooningCoefficients(SimpleNamespace):
    """Extended g, c, f and independent vector reconstruction diagnostics."""


@dataclass
class FitOptions:
    ntheta: int = 256
    dpsi: float = 0.03
    n_radial_fit: int = 9
    max_iterations: int = 200
    theta_tolerance: float = 1e-9
    coefficient_tolerance: float = 1e-10
    optimizer_tolerance: float = 1e-10
    n_phase_starts: int = 32
    n_x_starts: int = 33
    max_abs_delta: float = 0.98
    max_kappa: float = 5.0
    max_I0: float = 3.0
    verbose: bool = False


def _integer(value, name, minimum):
    if not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(value)


def _finite(*values):
    if not np.all(np.isfinite(values)):
        raise ValueError("All scalar parameters must be finite.")


def _vec(value, name="array", minimum=1):
    a = np.asarray(value, dtype=float)
    if a.ndim not in (1, 2) or (a.ndim == 2 and 1 not in a.shape):
        raise ValueError(f"{name} must be a vector.")
    a = a.reshape(-1)
    if a.size < minimum or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} needs at least {minimum} finite entries.")
    return a


def _closed_grid(theta):
    theta = _vec(theta, "theta", 4)
    if np.any(np.diff(theta) <= 0) or not np.isclose(
        theta[-1] - theta[0], _TWO_PI, rtol=0, atol=1e-8
    ):
        raise ValueError("theta must increase over one closed 2*pi period.")
    return theta


def _bounded(fun, lower, upper, tolerance):
    result = minimize_scalar(
        fun, bounds=(float(lower), float(upper)), method="bounded",
        options={"xatol": tolerance, "maxiter": 500},
    )
    if not result.success:
        raise RuntimeError(f"Bounded optimization failed: {result.message}")
    return float(result.x)


def _round_matlab(value):
    value = np.asarray(value)
    return np.sign(value) * np.floor(np.abs(value) + 0.5)


def _geometry(model, coeff, theta):
    if model == "cshape":
        A, B, C, G, H = coeff
        R = A + B * np.sin(theta) + C * np.cos(2 * theta)
        Z = G * np.cos(theta) - H * np.sin(2 * theta)
        Rt = B * np.cos(theta) - 2 * C * np.sin(2 * theta)
        Zt = -G * np.sin(theta) - 2 * H * np.cos(2 * theta)
        Rtt = -B * np.sin(theta) - 4 * C * np.cos(2 * theta)
        Ztt = -G * np.cos(theta) + 4 * H * np.sin(2 * theta)
    else:
        R0, a, x, b, d = coeff
        cs, sn = np.cos(theta), np.sin(theta)
        eta = theta + x * sn
        eta_t, eta_tt = 1 + x * cs, -x * sn
        h = ((1 + cs) / 2)**3
        ht = -(3 / 8) * (1 + cs)**2 * sn
        htt = (3 / 4) * (1 + cs) * sn**2 - (3 / 8) * (1 + cs)**2 * cs
        R, Z = R0 + a * np.cos(eta) - d * h, b * sn
        Rt, Zt = -a * np.sin(eta) * eta_t - d * ht, b * cs
        Rtt = -a * np.cos(eta) * eta_t**2 - a * np.sin(eta) * eta_tt - d * htt
        Ztt = -b * sn
    return R, Z, Rt, Zt, Rtt, Ztt


def _boundary(model, coeff, ntheta):
    theta = np.linspace(0, _TWO_PI, _integer(ntheta, "ntheta", 16))
    _finite(*coeff)
    return Boundary(theta, *_geometry(model, coeff, theta))


def circular_boundary(R0, r, *, ntheta=256):
    if r <= 0:
        raise ValueError("r must be positive.")
    return _boundary("indented", [R0, r, 0.0, r, 0.0], ntheta)


def miller_boundary(R0, r, delta, kappa, *, ntheta=256):
    _finite(R0, r, delta, kappa)
    if r <= 0 or kappa <= 0 or abs(delta) >= 1:
        raise ValueError("Require r > 0, kappa > 0, and |delta| < 1.")
    return _boundary(
        "indented", [R0, r, np.arcsin(delta), kappa * r, 0.0], ntheta
    )


def cshape_boundary(A, B, C, G, H, *, ntheta=256):
    return _boundary("cshape", [A, B, C, G, H], ntheta)


def indented_miller_boundary(R0, a, delta, kappa, I0, *, ntheta=256):
    _finite(R0, a, delta, kappa, I0)
    if a <= 0 or kappa <= 0 or abs(delta) >= 1 or I0 < 0:
        raise ValueError("Require a > 0, kappa > 0, |delta| < 1, and I0 >= 0.")
    bnd = _boundary(
        "indented", [R0, a, np.arcsin(delta), kappa * a, I0 * a], ntheta
    )
    bnd.x, bnd.b, bnd.d = np.arcsin(delta), kappa * a, I0 * a
    bnd.eta = bnd.theta + bnd.x * np.sin(bnd.theta)
    bnd.h = ((1 + np.cos(bnd.theta)) / 2)**3
    return bnd


def _prepare_contour(R, Z):
    R, Z = _vec(R, "R", 8), _vec(Z, "Z", 8)
    if R.size != Z.size:
        raise ValueError("R and Z must have equal lengths.")
    scale = max(1.0, np.max(np.hypot(R, Z)))
    if np.hypot(R[-1] - R[0], Z[-1] - Z[0]) <= 1e-12 * scale:
        R, Z = R[:-1], Z[:-1]
    if np.any(np.hypot(np.roll(R, -1) - R, np.roll(Z, -1) - Z)
              <= 16 * np.spacing(scale)):
        raise ValueError("The contour has coincident neighboring points.")
    return R.copy(), Z.copy()


def find_r_at_z0_from_surface(
    surface, *, z0=0.0, close_curve=True, tol_z=1e-12,
    unique_tol=1e-8, sort=True,
):
    R, Z = _vec(surface.R, "R", 2), _vec(surface.Z, "Z", 2)
    _finite(z0, tol_z, unique_tol)
    if R.size != Z.size or tol_z < 0 or unique_tol < 0:
        raise ValueError("Invalid contour sizes or tolerances.")
    if close_curve and np.hypot(R[-1] - R[0], Z[-1] - Z[0]) > 1e-14:
        R, Z = np.r_[R, R[0]], np.r_[Z, Z[0]]
    Z = Z - z0
    crossings = []
    for i in range(R.size - 1):
        if abs(Z[i]) <= tol_z:
            crossings.append(R[i])
        if Z[i] * Z[i + 1] < 0:
            crossings.append(R[i] - Z[i] * (R[i + 1] - R[i]) / (Z[i + 1] - Z[i]))
    if abs(Z[-1]) <= tol_z:
        crossings.append(R[-1])
    if not crossings:
        return np.empty(0)
    ordered = np.sort(crossings)
    groups = np.split(ordered, np.flatnonzero(np.diff(ordered) > unique_tol) + 1)
    merged = np.array([np.mean(group) for group in groups])
    if not sort:
        indices = [np.argmin(np.abs(np.asarray(crossings) - value)) for value in merged]
        merged = merged[np.argsort(indices)]
    return merged


def _arc_length_angle(R, Z):
    ds = np.hypot(np.diff(np.r_[R, R[0]]), np.diff(np.r_[Z, Z[0]]))
    if not np.isfinite(ds.sum()) or ds.sum() <= 0:
        raise ValueError("The contour has zero or invalid length.")
    return _TWO_PI * np.r_[0, np.cumsum(ds[:-1])] / ds.sum()


def _fit_circle(surface):
    R, Z = _vec(surface.R), _vec(surface.Z)
    crossings = find_r_at_z0_from_surface(surface)
    if crossings.size != 2:
        raise ValueError("The circle fit requires exactly two Z=0 crossings.")
    center, radius = crossings.mean(), np.ptp(crossings) / 2
    def objective(value):
        distance = np.hypot(R - value, Z)
        return np.mean((distance - distance.mean())**2)
    R0 = _bounded(objective, center - radius / 2, center + radius / 2, 1e-12)
    distances = np.hypot(R - R0, Z)
    r = distances.mean()
    return np.array([R0, r, np.sqrt(np.mean((distances - r)**2)),
                     np.max(np.abs(distances - r))])


def _fit_miller_surface(surface, ntheta):
    crossings = find_r_at_z0_from_surface(surface)
    if crossings.size != 2:
        raise ValueError("The Miller fit requires exactly two Z=0 crossings.")
    R0, r = crossings.mean(), np.ptp(crossings) / 2
    R, Z = _vec(surface.R), _vec(surface.Z)
    if np.hypot(R[-1] - R[0], Z[-1] - Z[0]) > 1e-12:
        R, Z = np.r_[R, R[0]], np.r_[Z, Z[0]]
    arc = np.r_[0, np.cumsum(np.hypot(np.diff(R), np.diff(Z)))]
    keep = np.r_[True, np.diff(arc) > 0]
    arc, R, Z = arc[keep], R[keep], Z[keep]
    rs, zs = CubicSpline(arc, R), CubicSpline(arc, Z)
    query = np.linspace(arc[0], arc[-1], max(4 * ntheta, 2048))
    i = int(np.argmax(zs(query)))
    top = query[i] if i in (0, query.size - 1) else _bounded(
        lambda x: -float(zs(x)), query[i - 1], query[i + 1], 1e-4
    )
    return np.array([R0, r, -(float(rs(top)) - R0) / r, float(zs(top)) / r])


def _project_to_ordered_angles(theta, separation):
    ramp = np.arange(theta.size) * separation
    levels, weights, starts, ends = [], [], [], []
    for i, value in enumerate(theta - ramp):
        levels.append(value); weights.append(1); starts.append(i); ends.append(i + 1)
        while len(levels) > 1 and levels[-2] > levels[-1]:
            total = weights[-2] + weights[-1]
            levels[-2] = (weights[-2] * levels[-2] + weights[-1] * levels[-1]) / total
            weights[-2], ends[-2] = total, ends[-1]
            levels.pop(); weights.pop(); starts.pop(); ends.pop()
    result = np.empty_like(theta)
    for value, start, end in zip(levels, starts, ends):
        result[start:end] = value
    return result + ramp


def _objective(model, coeff, theta, R, Z):
    rf, zf = _geometry(model, coeff, theta)[:2]
    return float(np.mean((R - rf)**2 + (Z - zf)**2))


def _closest_point_angles(model, coeff, theta_old, R, Z):
    indented = model == "indented"
    nsearch = max(256 if indented else 128,
                  (4 if indented else 2) * int(np.ceil(np.sqrt(R.size))))
    search = np.arange(nsearch) * (_TWO_PI / nsearch)
    rf, zf = _geometry(model, coeff, search)[:2]
    nearest = np.argmin((R[:, None] - rf)**2 + (Z[:, None] - zf)**2, axis=1)
    seed = search[nearest]
    theta = seed + _TWO_PI * _round_matlab((theta_old - seed) / _TWO_PI)
    max_step = 2 * _TWO_PI / nsearch
    for _ in range(15 if indented else 12):
        rf, zf, rt, zt, rtt, ztt = _geometry(model, coeff, theta)
        dr, dz = rf - R, zf - Z
        numerator = dr * rt + dz * zt
        denominator = rt**2 + zt**2 + dr * rtt + dz * ztt
        scale = rt**2 + zt**2 + np.abs(dr * rtt) + np.abs(dz * ztt)
        valid = np.abs(denominator) > 1e-12 * np.maximum(scale, _EPS)
        step = np.zeros(R.size)
        step[valid] = numerator[valid] / denominator[valid]
        step = np.clip(step, -max_step, max_step)
        best = dr**2 + dz**2
        step_scale, candidate = np.ones(R.size), theta.copy()
        for _ in range(8):
            trial = theta - step_scale * step
            rr, zz = _geometry(model, coeff, trial)[:2]
            distance = (rr - R)**2 + (zz - Z)**2
            improved = distance < best
            candidate[improved], best[improved] = trial[improved], distance[improved]
            step_scale[~improved] *= 0.5
        change = np.max(np.abs(candidate - theta))
        theta = candidate
        if change <= 1e-11:
            break
    separation = 64 * np.spacing(max(1.0, np.max(np.abs(theta))))
    ordered = _project_to_ordered_angles(theta, separation)
    previous = _objective(model, coeff, theta_old, R, Z)
    for backtrack in range(11):
        candidate = theta_old + 2.0**(-backtrack) * (ordered - theta_old)
        if np.all(np.diff(candidate) > 0) and candidate[-1] < candidate[0] + _TWO_PI:
            if _objective(model, coeff, candidate, R, Z) <= previous * (1 + 1e-12):
                return candidate
    return theta_old.copy()


def _linear_cshape(theta, R, Z):
    XR = np.column_stack((np.ones(theta.size), np.sin(theta), np.cos(2 * theta)))
    XZ = np.column_stack((np.cos(theta), -np.sin(2 * theta)))
    return np.r_[np.linalg.lstsq(XR, R, rcond=None)[0],
                 np.linalg.lstsq(XZ, Z, rcond=None)[0]]


def _valid_indented(coeff, options):
    R0, a, x, b, d = coeff
    return bool(
        np.all(np.isfinite(coeff)) and R0 > 0 and a > 0 and b > 0 and d >= 0
        and abs(np.sin(x)) <= options.max_abs_delta * (1 + 1e-12)
        and b / a <= options.max_kappa and d / a <= options.max_I0
    )


def _coefficients_at_x(x, theta, R, Z, options):
    h = ((1 + np.cos(theta)) / 2)**3
    eta = theta + x * np.sin(theta)
    XR = np.column_stack((np.ones(theta.size), np.cos(eta), -h))
    rc = np.linalg.lstsq(XR, R, rcond=None)[0]
    b = np.linalg.lstsq(np.sin(theta)[:, None], Z, rcond=None)[0][0]
    coeff = np.array([rc[0], rc[1], x, b, rc[2]])
    value = _objective("indented", coeff, theta, R, Z)
    if not _valid_indented(coeff, options):
        R0, a, _, b, d = coeff
        safe_a = max(abs(a), _EPS)
        penalty = sum(max(0.0, -v / safe_a)**2 for v in (R0, a, b, d))
        penalty += max(0.0, b / safe_a - options.max_kappa)**2
        penalty += max(0.0, d / safe_a - options.max_I0)**2
        scale = max(1.0, 0.25 * (np.ptp(R) + np.ptp(Z)))
        value += scale**2 * (100 + penalty)
    return coeff, value


def _fit_coefficients(model, theta, R, Z, options, x_seed=None):
    if model == "cshape":
        coeff = _linear_cshape(theta, R, Z)
        return coeff, _objective(model, coeff, theta, R, Z)
    bound = np.arcsin(options.max_abs_delta)
    grid = np.linspace(-bound, bound, options.n_x_starts)
    if x_seed is None:
        values = [_coefficients_at_x(x, theta, R, Z, options)[1] for x in grid]
        i = int(np.argmin(values))
        best = grid[i]
        lower, upper = grid[max(0, i - 1)], grid[min(grid.size - 1, i + 1)]
    else:
        best = np.clip(x_seed, -bound, bound)
        width = 4 * bound / (grid.size - 1)
        lower, upper = max(-bound, best - width), min(bound, best + width)
    coeff, value = _coefficients_at_x(best, theta, R, Z, options)
    if upper - lower > options.optimizer_tolerance:
        optimum = _bounded(
            lambda x: _coefficients_at_x(x, theta, R, Z, options)[1],
            lower, upper, options.optimizer_tolerance,
        )
        alt, alt_value = _coefficients_at_x(optimum, theta, R, Z, options)
        if alt_value < value:
            coeff, value = alt, alt_value
    return coeff, value


def _canonicalize_cshape(coeff):
    coeff = coeff.copy()
    B, G = coeff[1], coeff[3]
    if B >= 0 and G >= 0:
        transform = "identity"
    elif B < 0 and G < 0:
        coeff[[1, 3]] *= -1
        transform = "theta-plus-pi"
    elif B < 0:
        coeff[[1, 4]] *= -1
        transform = "theta-reflection"
    else:
        coeff[[3, 4]] *= -1
        transform = "pi-minus-theta"
    return coeff, transform


def _fit_nonlinear_surface(surface, model, options):
    R, Z = _prepare_contour(surface.R, surface.Z)
    best, best_objective = None, np.inf
    phase_step = _TWO_PI / options.n_phase_starts
    for orientation in (1, -1):
        rr, zz = (R, Z) if orientation == 1 else (R[::-1], Z[::-1])
        arc = _arc_length_angle(rr, zz)
        def phase_objective(phase):
            return _fit_coefficients(model, arc + phase, rr, zz, options)[1]
        phases = np.arange(options.n_phase_starts) * phase_step
        seed = phases[np.argmin([phase_objective(p) for p in phases])]
        phase = _bounded(
            phase_objective, seed - phase_step, seed + phase_step,
            options.optimizer_tolerance,
        )
        theta = arc + phase
        coeff, value = _fit_coefficients(model, theta, rr, zz, options)
        valid = model == "cshape" or _valid_indented(coeff, options)
        if valid and value < best_objective:
            best_objective = value
            best = (rr.copy(), zz.copy(), theta, coeff, orientation, phase)
    if best is None:
        raise ValueError("No admissible phase/orientation initialization was found.")
    R, Z, theta, coeff, orientation, phase = best
    initial_coeff = coeff.copy()
    history = [best_objective]
    converged = False
    theta_change = coefficient_change = np.inf
    for iteration in range(1, options.max_iterations + 1):
        updated = _closest_point_angles(model, coeff, theta, R, Z)
        cnew, value = _fit_coefficients(
            model, updated, R, Z, options,
            x_seed=coeff[2] if model == "indented" else None,
        )
        if model == "indented" and not _valid_indented(cnew, options):
            break
        history.append(value)
        theta_change = np.max(np.abs(updated - theta))
        coefficient_change = np.linalg.norm(cnew - coeff) / max(np.linalg.norm(coeff), _EPS)
        theta, coeff = updated, cnew
        if options.verbose:
            print(f"{model} {iteration}: RMS={np.sqrt(value):.6e}, "
                  f"dtheta={theta_change:.3e}, dcoeff={coefficient_change:.3e}")
        if theta_change <= options.theta_tolerance and coefficient_change <= options.coefficient_tolerance:
            converged = True
            break
    correspondence = _linear_cshape(theta, R, Z) if model == "cshape" else coeff
    rr, zz = _geometry(model, correspondence, theta)[:2]
    distance = np.hypot(rr - R, zz - Z)
    rms = np.sqrt(np.mean(distance**2))
    diagnostic = SimpleNamespace(
        converged=converged, iterations=iteration,
        rms_distance=float(rms), max_distance=float(distance.max()),
        relative_rms_distance=float(rms / max(0.25 * (np.ptp(R) + np.ptp(Z)), _EPS)),
        initial_rms_distance=float(np.sqrt(history[0])),
        objective_history=np.array(history),
        final_theta_change=float(theta_change),
        final_coefficient_change=float(coefficient_change),
        final_relative_objective_change=(
            abs(history[-1] - history[-2]) / max(history[-2], _EPS)
            if len(history) > 1 else np.nan
        ),
        theta=np.mod(theta, _TWO_PI), orientation=orientation,
        initial_phase=float(np.mod(phase, _TWO_PI)),
        initial_coeff=initial_coeff,
        vertical_center_offset=float(0.5 * (Z.max() + Z.min())),
    )
    if model == "cshape":
        coeff, diagnostic.canonical_transform = _canonicalize_cshape(coeff)
    if not converged and options.verbose:
        warnings.warn(f"{model} reached {iteration} iterations without convergence.",
                      RuntimeWarning, stacklevel=2)
    return coeff, diagnostic


def _resolve_input(eq, fields, psiN):
    if hasattr(eq, "data") and hasattr(eq, "fields"):
        if fields is not None and psiN is not None:
            raise TypeError("Use (eq, psiN) or (data, fields, psiN).")
        return eq.data, eq.fields, float(fields if psiN is None else psiN)
    if fields is None or psiN is None:
        raise TypeError("Use (eq, psiN) or (data, fields, psiN).")
    return eq, fields, float(psiN)


def _fit(eq, fields, psiN, model, options):
    data, fields, psiN = _resolve_input(eq, fields, psiN)
    _finite(psiN, options.dpsi)
    if not 0 < psiN < 1 or options.dpsi <= 0:
        raise ValueError("Require 0 < psiN < 1 and dpsi > 0.")
    _integer(options.ntheta, "ntheta", 16)
    _integer(options.n_radial_fit, "n_radial_fit", 3)
    if options.n_radial_fit % 2 == 0:
        raise ValueError("n_radial_fit must be odd.")
    if model in ("cshape", "indented"):
        _integer(options.max_iterations, "max_iterations", 1)
        _integer(options.n_phase_starts, "n_phase_starts", 8)
        tolerances = (options.theta_tolerance, options.coefficient_tolerance,
                      options.optimizer_tolerance)
        if not np.all(np.isfinite(tolerances)) or min(tolerances) <= 0:
            raise ValueError("Fitting tolerances must be positive and finite.")
    if model == "indented":
        _integer(options.n_x_starts, "n_x_starts", 9)
        _finite(options.max_abs_delta, options.max_kappa, options.max_I0)
        if not 0 < options.max_abs_delta < 1 or min(options.max_kappa, options.max_I0) <= 0:
            raise ValueError("Invalid indented-Miller bounds.")

    half_width = min(options.dpsi, 0.95 * min(psiN, 1 - psiN))
    labels = np.linspace(psiN - half_width, psiN + half_width, options.n_radial_fit)
    center = labels.size // 2
    central_surface = extract_flux_surface(data, fields, psiN, npoints=options.ntheta)
    coeffs, diagnostics = [], []
    for i, label in enumerate(labels):
        surf = central_surface if i == center else extract_flux_surface(
            data, fields, label, npoints=options.ntheta
        )
        if model == "circular":
            coeff, diagnostic = _fit_circle(surf), None
        elif model == "miller":
            coeff, diagnostic = _fit_miller_surface(surf, options.ntheta), None
        else:
            coeff, diagnostic = _fit_nonlinear_surface(surf, model, options)
        coeffs.append(coeff)
        diagnostics.append(diagnostic)
    coeffs = np.array(coeffs)
    coeff = coeffs[center]
    radius = coeff[1]
    offset = coeffs[:, 1] - radius
    if np.ptp(offset) <= 1e-10 * max(1.0, abs(radius)):
        raise ValueError("The fitted minor radius does not vary across surfaces.")
    design = np.column_stack((np.ones(labels.size), offset))
    regression_values = coeffs.copy()
    if model == "miller":
        regression_values[:, 2] = np.arcsin(np.clip(coeffs[:, 2], -1, 1))
    derivative = np.linalg.lstsq(design, regression_values, rcond=None)[0][1]
    physical_flux = data.simag + labels * (data.sibry - data.simag)
    dpsi_dr_fit = np.linalg.lstsq(design, physical_flux, rcond=None)[0][1]
    grid = np.linspace(0, 1, data.nw)
    F, q = float(np.interp(psiN, grid, data.fpol)), float(np.interp(psiN, grid, data.qpsi))
    radial = SimpleNamespace(
        psiN=labels, r=coeffs[:, 1], half_width=half_width,
        npoints=labels.size,
    )
    if model in ("circular", "miller"):
        R0, r = coeff[:2]
        radial.R0 = coeffs[:, 0]
        param = ShapeParameters(R0=R0, r=r, A=R0 / r, dR0_dr=derivative[0])
        if model == "circular":
            param.fit_residual_rms, param.fit_residual_max = coeff[2:]
            radial.residual_rms, radial.residual_max = coeffs[:, 2], coeffs[:, 3]
            bnd = circular_boundary(R0, r, ntheta=options.ntheta)
        else:
            delta, kappa = coeff[2:]
            param.delta, param.kappa = delta, kappa
            param.s_delta, param.s_kappa = r * derivative[2], r * derivative[3] / kappa
            radial.delta, radial.kappa = coeffs[:, 2], coeffs[:, 3]
            bnd = miller_boundary(R0, r, delta, kappa, ntheta=options.ntheta)
    elif model == "cshape":
        A, B, C, G, H = coeff
        R0 = A - C
        param = ShapeParameters(
            A=A, B=B, C=C, G=G, H=H, r=B, R0=R0, aspect_ratio=R0 / B,
            coeff=coeff.copy(), dcoeff_dr=derivative.copy(),
            coeff_table=np.column_stack((coeff, derivative)),
            fit=diagnostics[center], dpsi_dr_fit=dpsi_dr_fit,
        )
        for name, value in zip(("dA_dr", "dB_dr", "dC_dr", "dG_dr", "dH_dr"), derivative):
            setattr(param, name, value)
        radial.coeff, radial.dcoeff_dr, radial.diagnostics = coeffs, derivative, diagnostics
        bnd = cshape_boundary(*coeff, ntheta=options.ntheta)
    else:
        R0, a, x, b, d = coeff
        delta, kappa, I0 = np.sin(x), b / a, d / a
        param = ShapeParameters(
            R0=R0, a=a, r=a, A=R0 / a, x=x, delta=delta, b=b, kappa=kappa,
            d=d, I0=I0, indentation=I0, coeff=coeff.copy(),
            dR0_da=derivative[0], dR0_dr=derivative[0], da_da=derivative[1],
            dx_da=derivative[2], db_da=derivative[3], dd_da=derivative[4],
            dcoeff_da=derivative.copy(), coeff_table=np.column_stack((coeff, derivative)),
            s_delta=a * derivative[2], s_kappa=a * derivative[3] / b - 1,
            dI0_da=(derivative[4] - I0) / a,
            dI0_dr=(derivative[4] - I0) / a,
            s_I0=derivative[4] - I0, s_I=derivative[4] - I0,
            dpsi_da_fit=dpsi_dr_fit, dpsi_dr_fit=dpsi_dr_fit,
            fit=diagnostics[center],
        )
        radial.a, radial.R0, radial.x = coeffs[:, 1], coeffs[:, 0], coeffs[:, 2]
        radial.delta, radial.b, radial.d = np.sin(coeffs[:, 2]), coeffs[:, 3], coeffs[:, 4]
        radial.kappa, radial.I0 = coeffs[:, 3] / coeffs[:, 1], coeffs[:, 4] / coeffs[:, 1]
        radial.coeff, radial.dcoeff_da, radial.diagnostics = coeffs, derivative, diagnostics
        bnd = indented_miller_boundary(R0, a, delta, kappa, I0, ntheta=options.ntheta)
    param.F, param.q, param.B0, param.radial_fit = F, q, F / R0, radial
    bnd.Req, bnd.Zeq = central_surface.R, central_surface.Z
    return param, bnd


def fit_circular(eq, fields=None, psiN=None, *, ntheta=256, dpsi=0.03, n_radial_fit=9):
    return _fit(eq, fields, psiN, "circular", FitOptions(
        ntheta=ntheta, dpsi=dpsi, n_radial_fit=n_radial_fit
    ))


def fit_miller(eq, fields=None, psiN=None, *, ntheta=256, dpsi=0.03, n_radial_fit=9):
    return _fit(eq, fields, psiN, "miller", FitOptions(
        ntheta=ntheta, dpsi=dpsi, n_radial_fit=n_radial_fit
    ))


def fit_cshape(eq, fields=None, psiN=None, **options):
    return _fit(eq, fields, psiN, "cshape", FitOptions(**options))


def fit_indented_miller(eq, fields=None, psiN=None, **options):
    return _fit(eq, fields, psiN, "indented", FitOptions(**options))

# ----------------------------------------------------------------------
# Mercier-Luc geometry
# ----------------------------------------------------------------------

def _assemble_local(surface, bnd, Rr, Zr, Bp, dpsi_dr, *,
                    model, align_flux=False, jac=None):
    theta = _closed_grid(bnd.theta)
    R, Z, Rt, Zt, Rtt, Ztt = (
        _vec(getattr(bnd, name), name, 4)
        for name in ("R", "Z", "Rt", "Zt", "Rtt", "Ztt")
    )
    if any(v.size != theta.size for v in (R, Z, Rt, Zt, Rtt, Ztt, Rr, Zr, Bp)):
        raise ValueError("All local geometry arrays must have the same length.")
    dl = np.hypot(Rt, Zt)
    if np.any(R <= 0) or np.any(dl <= 0) or np.any(Bp <= 0):
        raise ValueError("R, tangent length, and Bp must be positive.")
    if not np.all(np.isfinite(np.r_[Rr, Zr, Bp])):
        raise ValueError("Local radial derivatives and Bp must be finite.")
    if not np.isfinite(surface.F) or not np.isfinite(surface.q) or surface.q == 0:
        raise ValueError("F and q must be finite and q must be nonzero.")
    if surface.F / surface.q <= 0:
        raise ValueError("F and q must have consistent signs for increasing PEST theta.")
    cosu_raw, sinu_raw = Rt / dl, -Zt / dl
    projection = sinu_raw * Rr + cosu_raw * Zr
    normal_sign = np.sign(np.median(projection))
    if align_flux:
        normal_sign *= np.sign(dpsi_dr)
    if normal_sign == 0:
        raise ValueError("The Mercier normal orientation is undefined.")
    raw_curvature = (Rt * Ztt - Zt * Rtt) / dl**3
    Bphi = surface.F / R
    B2 = Bp**2 + Bphi**2
    H = surface.F / (surface.q * R**2 * Bp)
    theta_PEST = cumulative_trapezoid(H * dl, theta, initial=0)
    return LocalEquilibrium(
        model=model, theta_geo=theta, theta_PEST=theta_PEST,
        R=R, Z=Z, Rt=Rt, Zt=Zt, Rtt=Rtt, Ztt=Ztt,
        Rr=np.asarray(Rr), Zr=np.asarray(Zr),
        jac=Rr * Zt - Rt * Zr if jac is None else jac,
        dl_dt=dl, cosu=normal_sign * cosu_raw, sinu=normal_sign * sinu_raw,
        cosu_raw=cosu_raw, sinu_raw=sinu_raw, normal_sign=float(normal_sign),
        normal_projection_raw=projection,
        invRc=normal_sign * raw_curvature, invRc_raw=raw_curvature,
        dpsi_dr=float(dpsi_dr), Bp=Bp, Bphi=Bphi, B2=B2, B=np.sqrt(B2),
        H=H, F=float(surface.F), q=float(surface.q),
        q_check=float(surface.F / _TWO_PI * trapezoid(dl / (R**2 * Bp), theta)),
    )


def circular_mercier_luc(surface, param, bnd):
    theta = _closed_grid(bnd.theta)
    factor = 1 + param.dR0_dr * np.cos(theta)
    if np.any(factor <= 0) or param.R0 <= param.r:
        raise ValueError("Require R0 > r and 1 + dR0_dr*cos(theta) > 0.")
    Rr, Zr = param.dR0_dr + np.cos(theta), np.sin(theta)
    jac = param.r * factor
    D = np.sqrt(param.R0**2 - param.r**2)
    exact = _TWO_PI * (param.r / D + param.dR0_dr * (1 - param.R0 / D))
    numeric = trapezoid(jac / bnd.R, theta)
    dpsi_dr = surface.F * exact / (_TWO_PI * surface.q)
    Bp = abs(dpsi_dr) / (bnd.R * factor)
    result = _assemble_local(
        surface, bnd, Rr, Zr, Bp, dpsi_dr, model="circular", jac=jac
    )
    result.radial_metric_factor = factor
    result.q_integral_analytic, result.q_integral_numeric = float(exact), float(numeric)
    result.q_integral_relative_error = float(abs(numeric - exact) / max(abs(exact), _EPS))
    actual_jac = Rr * bnd.Zt - bnd.Rt * Zr
    result.jac_relative_error = float(np.max(np.abs(actual_jac - jac)) / max(np.max(np.abs(jac)), _EPS))
    return result


def _analytic_poloidal_field(surface, bnd, Rr, Zr, flux_sign=None):
    jac = Rr * bnd.Zt - bnd.Rt * Zr
    scale = max(1.0, np.max(np.abs(Rr * bnd.Zt)) + np.max(np.abs(bnd.Rt * Zr)))
    if not np.all(np.isfinite(jac)) or np.any(np.abs(jac) <= 1e-12 * scale):
        raise ValueError("The fitted radial coordinate has a singular Jacobian.")
    if surface.q == 0 or not np.isfinite(surface.q):
        raise ValueError("Surface q must be finite and nonzero.")
    integral = trapezoid(np.abs(jac) / bnd.R, bnd.theta)
    dpsi = surface.F * integral / (_TWO_PI * surface.q)
    if flux_sign is not None:
        if not np.isfinite(flux_sign) or flux_sign == 0:
            raise ValueError("The fitted flux direction is undefined.")
        dpsi = np.sign(flux_sign) * abs(dpsi)
    Bp = abs(dpsi) * np.hypot(bnd.Rt, bnd.Zt) / (np.abs(jac) * bnd.R)
    return Bp, dpsi, jac


def miller_mercier_luc(surface, param, bnd):
    theta, r = bnd.theta, param.r
    eta = theta + np.arcsin(param.delta) * np.sin(theta)
    Rr = param.dR0_dr + np.cos(eta) - np.sin(eta) * param.s_delta * np.sin(theta)
    Zr = param.kappa * (1 + param.s_kappa) * np.sin(theta)
    Bp, dpsi, jac = _analytic_poloidal_field(surface, bnd, Rr, Zr)
    return _assemble_local(
        surface, bnd, Rr, Zr, Bp, dpsi, model="miller", jac=jac
    )


def cshape_mercier_luc(surface, param, bnd):
    theta = bnd.theta
    Rr = param.dA_dr + param.dB_dr * np.sin(theta) + param.dC_dr * np.cos(2 * theta)
    Zr = param.dG_dr * np.cos(theta) - param.dH_dr * np.sin(2 * theta)
    Bp, dpsi, jac = _analytic_poloidal_field(
        surface, bnd, Rr, Zr, flux_sign=param.dpsi_dr_fit
    )
    return _assemble_local(
        surface, bnd, Rr, Zr, Bp, dpsi, model="cshape", align_flux=True, jac=jac
    )


def indented_miller_mercier_luc(surface, param, bnd):
    theta = bnd.theta
    eta = theta + param.x * np.sin(theta)
    h = ((1 + np.cos(theta)) / 2)**3
    Ra = (
        param.dR0_da + param.da_da * np.cos(eta)
        - param.a * np.sin(eta) * param.dx_da * np.sin(theta) - param.dd_da * h
    )
    Za = param.db_da * np.sin(theta)
    Bp, dpsi, jac = _analytic_poloidal_field(
        surface, bnd, Ra, Za, flux_sign=param.dpsi_da_fit
    )
    result = _assemble_local(
        surface, bnd, Ra, Za, Bp, dpsi,
        model="indented-miller", align_flux=True, jac=jac,
    )
    result.Ra, result.Za, result.dpsi_da = Ra, Za, dpsi
    return result


def _fourier_geometry(value, maximum_mode):
    value = _vec(value)
    n = value.size
    modes = np.r_[np.arange((n + 1) // 2), -np.arange(n // 2, 0, -1)]
    if n % 2 == 0:
        # Match the original even-grid Nyquist convention.
        modes[n // 2] = 0
    spectrum = np.fft.fft(value)
    spectrum[np.abs(modes) > maximum_mode] = 0
    return (
        np.fft.ifft(spectrum).real,
        np.fft.ifft(1j * modes * spectrum).real,
        np.fft.ifft(-modes**2 * spectrum).real,
    )


def _close(core):
    return np.r_[core, core[0]]


def flux_surface_mercier_luc(
    surface, fields, *, fourier_modes=32, projection_iterations=2, enforce_q=True,
):
    R, Z = _prepare_contour(surface.R, surface.Z)
    n = R.size
    if n < 32:
        raise ValueError("The direct model requires at least 32 distinct points.")
    _integer(fourier_modes, "fourier_modes", 2)
    _integer(projection_iterations, "projection_iterations", 0)
    if fourier_modes > (n - 1) // 2 or projection_iterations > 10:
        raise ValueError("Too many Fourier modes or projection iterations.")
    area = np.sum(R * np.roll(Z, -1) - np.roll(R, -1) * Z) / 2
    if area <= 0:
        raise ValueError("The direct model requires a counter-clockwise contour.")
    for _ in range(projection_iterations):
        Rf = _fourier_geometry(R, fourier_modes)[0]
        Zf = _fourier_geometry(Z, fourier_modes)[0]
        value = fields.eval(Rf, Zf)
        grad2 = value.psiR**2 + value.psiZ**2
        if not np.all(np.isfinite(grad2)) or np.any(grad2 <= _EPS):
            raise ValueError("Invalid grad(psi) on the projected contour.")
        correction = (surface.psi - value.psi) / grad2
        R, Z = Rf + correction * value.psiR, Zf + correction * value.psiZ
    R, Rt, Rtt = _fourier_geometry(R, fourier_modes)
    Z, Zt, Ztt = _fourier_geometry(Z, fourier_modes)
    value = fields.eval(R, Z)
    grad = _vec(value.gradPsi)
    dl = np.hypot(Rt, Zt)
    if np.any(grad <= 0) or np.any(dl <= 0) or np.any(R <= 0):
        raise ValueError("The direct geometry contains an invalid gradient or tangent.")
    raw_projection = (-Zt * value.psiR + Rt * value.psiZ) / dl
    normal_sign = np.sign(np.median(raw_projection))
    aligned = normal_sign * raw_projection
    if normal_sign == 0 or np.any(aligned <= 0):
        raise ValueError("The smooth contour normal is not aligned with grad(psi).")
    Bp_raw = grad / R
    Bp = np.exp(_fourier_geometry(np.log(Bp_raw), fourier_modes)[0])
    theta = np.linspace(0, _TWO_PI, n + 1)
    bnd = Boundary(theta, *(_close(v) for v in (R, Z, Rt, Zt, Rtt, Ztt)))
    q_before = surface.F / _TWO_PI * trapezoid(
        _close(dl) / (bnd.R**2 * _close(Bp)), theta
    )
    correction = q_before / surface.q
    if not np.isfinite(correction) or correction <= 0:
        raise ValueError("The signs of F and q do not define a positive Bp correction.")
    if enforce_q:
        Bp = correction * Bp
    else:
        correction = 1.0
    grad2 = value.psiR**2 + value.psiZ**2
    Rpsi, Zpsi = _close(value.psiR / grad2), _close(value.psiZ / grad2)
    result = _assemble_local(
        surface, bnd, Rpsi, Zpsi, _close(Bp), 1.0, model="direct-flux-surface"
    )
    result.Rpsi, result.Zpsi = Rpsi, Zpsi
    result.normal_projection_raw = _close(raw_projection)
    result.diagnostics = SimpleNamespace(
        FourierModes=fourier_modes, ProjectionIterations=projection_iterations,
        qBeforeCorrection=float(q_before), qCorrection=float(correction),
        qRelativeError=float((result.q_check - surface.q) / max(abs(surface.q), _EPS)),
        psiNRms=float(np.sqrt(np.mean((value.psiN - surface.psiN)**2))),
        psiNMaximum=float(np.max(np.abs(value.psiN - surface.psiN))),
        minimumNormalAlignment=float(np.min(aligned / grad)),
        PESTClosure=float(result.theta_PEST[-1] - _TWO_PI),
        BpRaw=_close(Bp_raw),
    )
    return result


# ----------------------------------------------------------------------
# Radial derivative of PEST theta
# ----------------------------------------------------------------------

def _periodic_derivative_closed(x, value):
    x = _closed_grid(x)
    value = _vec(value)
    if value.size != x.size or not np.allclose(
        np.diff(x), (x[-1] - x[0]) / (x.size - 1), rtol=1e-10, atol=1e-12
    ):
        raise ValueError("This periodic derivative requires a uniform closed grid.")
    step = (x[-1] - x[0]) / (x.size - 1)
    core = value[:-1]
    return _close((np.roll(core, -1) - np.roll(core, 1)) / (2 * step))


def evaluate_theta_psi(merluc, FFprime, pprime, qprime, *, enforce_periodicity=True):
    _finite(FFprime, pprime, qprime)
    theta = _closed_grid(merluc.theta_geo)
    pest = _vec(merluc.theta_PEST, "theta_PEST", 4)
    if pest.size != theta.size or np.any(np.diff(pest) <= 0):
        raise ValueError("theta_PEST must be increasing on the same closed grid.")
    R, Bp, H, dl = merluc.R, merluc.Bp, merluc.H, merluc.dl_dt
    F, q = merluc.F, merluc.q
    if F == 0 or q == 0 or np.any(Bp <= 0):
        raise ValueError("F and q must be nonzero and Bp must be positive.")
    base = -qprime / q + (
        -2 * merluc.invRc - 2 * merluc.sinu / R + _MU0 * pprime * R / Bp
    ) / (R * Bp)
    coefficient = 1 / F**2 + 1 / (R**2 * Bp**2)
    weight = H * dl
    denominator = trapezoid(weight * coefficient, theta)
    if abs(denominator) <= _EPS:
        raise ValueError("The FFprime periodicity coefficient is numerically zero.")
    FFprime_periodic = -trapezoid(weight * base, theta) / denominator
    used = FFprime_periodic if enforce_periodicity else FFprime
    Gamma_input, Gamma = base + FFprime * coefficient, base + used * coefficient
    input_derivative = cumulative_trapezoid(weight * Gamma_input, theta, initial=0)
    theta_psi = cumulative_trapezoid(weight * Gamma, theta, initial=0)
    pest_derivative = cumulative_trapezoid(Gamma, pest, initial=0)
    diagnostic = ThetaDiagnostic(
        model=merluc.model, theta=pest, theta_geo=theta,
        theta_psi=theta_psi, Gamma=Gamma,
        R=R, Z=merluc.Z, Bp=Bp, Bphi=merluc.Bphi,
        B2=merluc.B2, B=merluc.B, H=H, dl_dt=dl,
        normal_R=merluc.sinu, normal_Z=merluc.cosu,
        tangent_R=merluc.Rt / dl, tangent_Z=merluc.Zt / dl,
        curvature_normal=merluc.invRc,
        dBp_dl=_periodic_derivative_closed(theta, Bp) / dl,
        F=F, q=q, normal_sign=merluc.normal_sign,
        FFprime_input=float(FFprime), FFprime_periodic=float(FFprime_periodic),
        FFprime_used=float(used), Gamma_input=Gamma_input,
        theta_psi_input=input_derivative, theta_psi_PEST=pest_derivative,
        closure_input=float(input_derivative[-1] - input_derivative[0]),
        closure_used=float(theta_psi[-1] - theta_psi[0]),
        closure_PEST=float(pest_derivative[-1] - pest_derivative[0]),
        integration_max_abs_difference=float(np.max(np.abs(theta_psi - pest_derivative))),
        input_closure_defect=float(
            abs(input_derivative[-1] - input_derivative[0])
            / max(trapezoid(np.abs(weight * Gamma_input), theta), _EPS)
        ),
        periodicity_enforced=bool(enforce_periodicity),
        FFprime=float(used), pprime=float(pprime), qprime=float(qprime),
    )
    return Gamma, theta_psi, diagnostic

# ----------------------------------------------------------------------
# Extended ballooning coefficients
# ----------------------------------------------------------------------

def _extend_periodic(theta, value, query):
    value = np.asarray(value, dtype=float)
    if value.ndim != 1 or value.size != theta.size or not np.all(np.isfinite(value)):
        raise ValueError("A periodic field must have one finite value per theta.")
    period = theta[-1] - theta[0]
    core = theta[:-1]
    grid = np.r_[core - period, core, core + period]
    data = np.tile(value[:-1], 3)
    return PchipInterpolator(grid, data)(np.mod(query - theta[0], period) + theta[0])


def _relative_error(actual, expected):
    return float(np.max(np.abs(actual - expected)) / max(np.max(np.abs(expected)), _EPS))


def evaluate_gcf(diag, theta_bal, *, theta_bnd=5 * np.pi, B_N=1.0, a_N=1.0):
    """Extend periodic geometry while retaining the secular qprime*theta term."""
    _finite(theta_bal, theta_bnd, B_N, a_N)
    if min(theta_bnd, B_N, a_N) <= 0:
        raise ValueError("theta_bnd, B_N, and a_N must be positive.")
    theta = _closed_grid(diag.theta)
    step = _TWO_PI / (theta.size - 1)
    nhalf = max(1, int(np.floor(theta_bnd / step + 0.5)))
    extended = np.arange(-nhalf, nhalf + 1) * step
    names = (
        "R", "Z", "Bp", "B2", "H", "theta_psi", "normal_R", "normal_Z",
        "tangent_R", "tangent_Z", "curvature_normal", "dBp_dl",
    )
    values = {name: _extend_periodic(theta, getattr(diag, name), extended) for name in names}
    R, Z, Bp = (values[name] for name in ("R", "Z", "Bp"))
    F, q, qprime, pprime = diag.F, diag.q, diag.qprime, diag.pprime
    if F == 0 or q == 0 or np.any(R <= 0) or np.any(Bp <= 0):
        raise ValueError("Invalid surface constants or extended magnetic fields.")
    Bphi = F / R
    B2 = Bp**2 + Bphi**2
    B = np.sqrt(B2)
    H = F / (q * R**2 * Bp)
    normal_norm = np.hypot(values["normal_R"], values["normal_Z"])
    tangent_norm = np.hypot(values["tangent_R"], values["tangent_Z"])
    if np.any(normal_norm == 0) or np.any(tangent_norm == 0):
        raise ValueError("The interpolated normal or tangent is degenerate.")
    nr, nz = values["normal_R"] / normal_norm, values["normal_Z"] / normal_norm
    tr, tz = values["tangent_R"] / tangent_norm, values["tangent_Z"] / tangent_norm
    theta_psi = values["theta_psi"]
    curvature, dBp_dl = values["curvature_normal"], values["dBp_dl"]
    delta = extended - theta_bal
    I = qprime * delta + q * theta_psi
    K = 1 / R**2 + F**2 / (R**4 * Bp**2) + (R * Bp)**2 * I**2
    P_normal = Bp**2 * curvature - F**2 * nr / R**3
    P_tangent = Bp * dBp_dl - F**2 * tr / R**3
    D = B2 / (R * Bp) * P_normal - F * Bp * I * P_tangent

    dpbar_dpsiN = _MU0 * a_N**2 * pprime / B_N
    dq_dpsiN = a_N**2 * B_N * qprime
    BoverBN = B / B_N
    b_dot_gradN_theta = a_N * F / (q * R**2 * B)
    gradN_alpha_sq = a_N**2 * K
    g = b_dot_gradN_theta * gradN_alpha_sq / BoverBN
    f = gradN_alpha_sq / (b_dot_gradN_theta * BoverBN**3)
    c = 2 * a_N * q * B_N**2 * R**2 * dpbar_dpsiN * D / (F * B2**2)

    # Independent reconstruction in the orthonormal [e_R, e_phi, e_Z] basis.
    zeros = np.zeros_like(R)
    normal = np.column_stack((nr, zeros, nz))
    tangent = np.column_stack((tr, zeros, tz))
    grad_psi = (R * Bp)[:, None] * normal
    grad_theta = (R * Bp * theta_psi)[:, None] * normal + H[:, None] * tangent
    grad_zeta = np.column_stack((zeros, 1 / R, zeros))
    grad_alpha = grad_zeta - q * grad_theta - (qprime * delta)[:, None] * grad_psi
    Bvec = np.column_stack((Bp * tr, Bphi, Bp * tz))
    bvec = Bvec / B[:, None]
    grad_P = P_normal[:, None] * normal + P_tangent[:, None] * tangent
    D_vector = np.sum(np.cross(Bvec, grad_P) * grad_alpha, axis=1)
    gradN_alpha = a_N * grad_alpha
    gradN_H = (a_N / B_N**2) * grad_P
    gradN_alpha_sq_vector = np.sum(gradN_alpha**2, axis=1)
    b_dot_gradN_theta_vector = np.sum(bvec * (a_N * grad_theta), axis=1)
    cross_term = np.sum(np.cross(bvec, gradN_H) * gradN_alpha, axis=1)
    g_vector = b_dot_gradN_theta_vector * gradN_alpha_sq_vector / BoverBN
    f_vector = gradN_alpha_sq_vector / (b_dot_gradN_theta_vector * BoverBN**3)
    c_vector = 2 / b_dot_gradN_theta_vector / BoverBN**4 * dpbar_dpsiN * cross_term
    validation = SimpleNamespace(
        basis_normal_norm_error=float(np.max(np.abs(np.sum(normal**2, axis=1) - 1))),
        basis_tangent_norm_error=float(np.max(np.abs(np.sum(tangent**2, axis=1) - 1))),
        basis_orthogonality_error=float(np.max(np.abs(np.sum(normal * tangent, axis=1)))),
        B2_interpolation_relative_error=_relative_error(values["B2"], B2),
        H_interpolation_relative_error=_relative_error(values["H"], H),
        K_relative_error=_relative_error(K, np.sum(grad_alpha**2, axis=1)),
        D_relative_error=_relative_error(D, D_vector),
        b_dot_gradN_theta_relative_error=_relative_error(
            b_dot_gradN_theta, b_dot_gradN_theta_vector
        ),
        g_relative_error=_relative_error(g, g_vector),
        f_relative_error=_relative_error(f, f_vector),
        c_relative_error=_relative_error(c, c_vector),
        all_finite=bool(np.all(np.isfinite(np.r_[g, c, f, I, K, D]))),
    )
    validation.vector_consistent = bool(max(
        validation.K_relative_error, validation.D_relative_error,
        validation.b_dot_gradN_theta_relative_error,
        validation.g_relative_error, validation.f_relative_error,
        validation.c_relative_error,
    ) <= 1e-8)
    if not validation.vector_consistent:
        warnings.warn(
            "Scalar/vector ballooning formulas disagree for this normal orientation. "
            "MATLAB scalar coefficients are retained; inspect validation before use.",
            RuntimeWarning, stacklevel=2,
        )
    return BallooningCoefficients(
        model=getattr(diag, "model", "miller-mercier-luc"),
        theta=extended, theta0=float(theta_bal),
        theta_bnd_requested=float(theta_bnd), theta_bnd_actual=float(np.max(np.abs(extended))),
        g=g, c=c, f=f, I=I, K=K, D=D, P_normal=P_normal, P_tangent=P_tangent,
        R=R, Z=Z, Bp=Bp, Bphi=Bphi, B2=B2, B=B, BoverBN=BoverBN, H=H,
        theta_psi=theta_psi, normal=normal, tangent=tangent,
        curvature_normal=curvature, dBp_dl=dBp_dl,
        F=F, q=q, qprime=qprime, pprime=pprime, FFprime=diag.FFprime,
        aN=a_N, BN=B_N, dq_dpsiN=dq_dpsiN, dpbar_dpsiN=dpbar_dpsiN,
        grad_psi=grad_psi, grad_theta=grad_theta, grad_zeta=grad_zeta,
        grad_alpha=grad_alpha, gradN_alpha=gradN_alpha, gradN_alpha_sq=gradN_alpha_sq,
        b_dot_gradN_theta=b_dot_gradN_theta, Bvec=Bvec, bvec=bvec,
        grad_P=grad_P, gradN_H=gradN_H, cross_term=cross_term,
        D_vector=D_vector, g_vector=g_vector, c_vector=c_vector, f_vector=f_vector,
        validation=validation,
    )


if __name__ == "__main__":
    import argparse
    from eq_process import Equilibrium

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geqdsk")
    parser.add_argument("--model", choices=("circular", "miller", "cshape", "indented", "direct"),
                        default="miller")
    parser.add_argument("--psiN", type=float, default=0.7)
    parser.add_argument("--ntheta", type=int, default=256)
    parser.add_argument("--max-iterations", type=int, default=200)
    args = parser.parse_args()
    eq = Equilibrium.from_geqdsk(args.geqdsk)
    surf = eq.extract_flux_surface(args.psiN, npoints=max(1024, args.ntheta))
    if args.model == "direct":
        local = flux_surface_mercier_luc(surf, eq.fields)
    else:
        fitters = dict(circular=fit_circular, miller=fit_miller,
                       cshape=fit_cshape, indented=fit_indented_miller)
        builders = dict(circular=circular_mercier_luc, miller=miller_mercier_luc,
                        cshape=cshape_mercier_luc, indented=indented_miller_mercier_luc)
        options = {"ntheta": args.ntheta}
        if args.model in ("cshape", "indented"):
            options["max_iterations"] = args.max_iterations
        param, bnd = fitters[args.model](eq, args.psiN, **options)
        local = builders[args.model](surf, param, bnd)
        if hasattr(param, "fit"):
            print(f"Fit converged: {param.fit.converged}; "
                  f"relative RMS: {param.fit.relative_rms_distance:.6e}")
    Gamma, theta_psi, diagnostic = evaluate_theta_psi(
        local, float(eq.fields.FFprime(args.psiN)),
        float(eq.fields.pprime(args.psiN)), float(eq.fields.qprime(args.psiN)),
    )
    coeff = evaluate_gcf(diagnostic, 0.0, a_N=eq.data.rmaxis, B_N=abs(eq.data.bcentr))
    print(f"q target/check: {local.q:.12g} / {local.q_check:.12g}")
    print(f"theta_psi closure: {diagnostic.closure_used:.6e}")
    print(f"g/c/f finite: {coeff.validation.all_finite}")
    print(f"Vector reconstruction errors (g,c,f): "
          f"{coeff.validation.g_relative_error:.3e}, "
          f"{coeff.validation.c_relative_error:.3e}, "
          f"{coeff.validation.f_relative_error:.3e}")
