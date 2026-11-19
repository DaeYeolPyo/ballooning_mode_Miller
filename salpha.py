"""Local-equilibrium s-alpha preparation, scans, marginal roots, and plots.

Dependencies: eq_process.py, local_equilibrium.py, sturm_liouville.py.
Positive lambdaMax denotes instability. Negative eigenvalues are retained.
Public function options use snake_case; result fields follow the MATLAB names.
"""

from dataclasses import replace
from types import SimpleNamespace
from pathlib import Path
import hashlib
import inspect
import json
import os
import time
import uuid
import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import trapezoid

from eq_process import Equilibrium, FluxSurface, SurfaceIntegrals
from local_equilibrium import (
    Boundary,
    fit_circular, fit_miller, fit_cshape, fit_indented_miller,
    circular_mercier_luc, miller_mercier_luc,
    cshape_mercier_luc, indented_miller_mercier_luc,
    flux_surface_mercier_luc, evaluate_theta_psi, evaluate_gcf,
)
from sturm_liouville import Sturm_Liouville

_MU0 = 4 * np.pi * 1e-7
_TWO_PI = 2 * np.pi
_CHECKPOINT_VERSION = 1


class SAlphaState(SimpleNamespace):
    """Fixed surface geometry and equilibrium-dependent preprocessing."""


class SAlphaPoint(SimpleNamespace):
    """One (shat, alpha), retaining every requested ballooning phase."""


class SAlphaScanResult(SimpleNamespace):
    """Scan arrays use [shat, alpha, theta0] axis order."""


class RefinedSAlphaResult(SAlphaScanResult):
    """Discovery scan plus freshly evaluated marginal-stability roots."""


class MarginalRoot(SimpleNamespace):
    """One transition with its final alpha and eigenvalue brackets."""


class MarginalTrace(SimpleNamespace):
    """One stable-to-unstable bracket per requested shat."""


def _vector(value, name, minimum=1):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(f"{name} must be real.")
    a = np.asarray(raw, dtype=float)
    if a.ndim == 0:
        a = a.reshape(1)
    if a.ndim not in (1, 2) or (a.ndim == 2 and 1 not in a.shape):
        raise ValueError(f"{name} must be a vector.")
    a = a.reshape(-1)
    if a.size < minimum or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} needs at least {minimum} finite entries.")
    return a.copy()


def _axis(value, name, *, minimum=2, nonnegative=False):
    a = _vector(value, name, minimum)
    if np.any(np.diff(a) <= 0) or (nonnegative and np.any(a < 0)):
        raise ValueError(f"{name} must be strictly increasing"
                         + (" and nonnegative." if nonnegative else "."))
    return a


def _integer(value, name, minimum):
    if not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(value)


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite.")
    return value


def _point_options(theta0, theta_boundary, n_balloon_points,
                   enforce_periodicity, require_vector_consistency):
    if theta0 is None or np.asarray(theta0).size == 0:
        theta0 = _TWO_PI * (np.arange(12) - 6) / 12
    theta0 = _vector(theta0, "theta0")
    if np.unique(theta0).size != theta0.size:
        raise ValueError("theta0 must contain distinct values.")
    boundary = _positive(theta_boundary, "theta_boundary")
    count = _integer(n_balloon_points, "n_balloon_points", 21)
    if count % 2 == 0:
        raise ValueError("n_balloon_points must be odd so theta=0 is represented.")
    return SimpleNamespace(
        theta0=theta0, thetaBoundary=boundary, nBalloonPoints=count,
        thetaGrid=np.linspace(-boundary, boundary, count),
        enforcePeriodicity=bool(enforce_periodicity),
        requireVectorConsistency=bool(require_vector_consistency),
    )


def _surface_constants(state):
    V, Vprime = float(state.integrals.V), float(state.integrals.Vprime)
    q, Raxis = float(state.surface.q), float(state.surface.Raxis)
    if (
        not np.all(np.isfinite([V, Vprime, q, Raxis]))
        or V <= 0 or Vprime == 0 or q == 0 or Raxis <= 0
    ):
        raise ValueError("Require V > 0, signed Vprime != 0, q != 0, Raxis > 0.")
    _positive(state.normalization.aN, "normalization.aN")
    _positive(state.normalization.BN, "normalization.BN")
    return V, Vprime, q, Raxis


def shat_to_qprime(shat, state):
    V, Vprime, q, _ = _surface_constants(state)
    return np.asarray(shat, dtype=float) * q * Vprime / (2 * V)


def alpha_to_pprime(alpha, state):
    V, Vprime, _, Raxis = _surface_constants(state)
    return (-np.asarray(alpha, dtype=float) * (4 * np.pi**2 / (2 * Vprime))
            * np.sqrt(2 * np.pi**2 * Raxis / V) / _MU0)


def local_volume_integrals(local, parameters=None):
    """Return V, signed Vprime and dV/dr (r=a for indented Miller)."""
    theta, R, Zt = local.theta_geo, local.R, local.Zt
    V = float(np.pi * abs(trapezoid(R**2 * Zt, theta)))
    dVdr = float(_TWO_PI * trapezoid(R * np.abs(local.jac), theta))
    if local.model == "circular" and parameters is not None:
        R0, r, shift = parameters.R0, parameters.r, parameters.dR0_dr
        V = float(2 * np.pi**2 * R0 * r**2)
        dVdr = float(2 * np.pi**2 * (2 * R0 * r + shift * r**2))
    derivative = float(local.dpsi_dr)
    if not np.isfinite(derivative) or abs(derivative) <= np.finfo(float).eps:
        raise ValueError("The signed local flux derivative is zero or invalid.")
    Vprime = dVdr / derivative
    if not np.all(np.isfinite([V, Vprime, dVdr])) or min(V, dVdr) <= 0 or Vprime == 0:
        raise ValueError("Invalid local volume integrals.")
    return SimpleNamespace(V=V, Vprime=float(Vprime), dVdr=dVdr)


def _resolve_equilibrium(eq, fields, target):
    if isinstance(eq, (str, Path)):
        eq = Equilibrium.from_geqdsk(eq)
    if hasattr(eq, "data") and hasattr(eq, "fields"):
        if fields is not None and target is not None:
            raise TypeError("Use (eq, target_psiN) or (data, fields, target_psiN).")
        target = fields if target is None else target
        return eq, float(target)
    if fields is None or target is None:
        raise TypeError("Use (eq, target_psiN) or (data, fields, target_psiN).")
    return Equilibrium(eq, fields), float(target)


def prepare_salpha_equilibrium(
    eq, fields=None, target_psiN=None, *, model="miller",
    n_surface_points=1024, n_fit_points=None,
    radial_fit_half_width=0.03, n_radial_fit=9,
    enforce_periodicity=True, fourier_modes=32,
    projection_iterations=2, enforce_direct_q=True,
    fit_options=None, volume_source="auto",
):
    """Prepare fixed geometry; fit_options holds nonlinear fitter options."""
    eq, target = _resolve_equilibrium(eq, fields, target_psiN)
    if not np.isfinite(target) or not 0 < target < 1:
        raise ValueError("target_psiN must be strictly between 0 and 1.")
    _integer(n_surface_points, "n_surface_points", 32)
    fit_points_given = n_fit_points is not None
    n_fit_points = n_surface_points if n_fit_points is None else n_fit_points
    _integer(n_fit_points, "n_fit_points", 16)
    width = _positive(radial_fit_half_width, "radial_fit_half_width")
    _integer(n_radial_fit, "n_radial_fit", 3)
    if n_radial_fit % 2 == 0:
        raise ValueError("n_radial_fit must be odd.")
    name = str(model).strip().lower().replace("_", "-")
    name = {"circle": "circular", "c-shape": "cshape", "flux-surface": "direct",
            "indented": "indented-miller"}.get(name, name)
    fitters = {
        "circular": (fit_circular, circular_mercier_luc),
        "miller": (fit_miller, miller_mercier_luc),
        "cshape": (fit_cshape, cshape_mercier_luc),
        "indented-miller": (fit_indented_miller, indented_miller_mercier_luc),
    }
    if name not in (*fitters, "direct"):
        raise ValueError("model must be circular, miller, cshape, indented-miller, or direct.")
    if not fit_points_given and name in ("cshape", "indented-miller"):
        n_fit_points = min(n_surface_points, 256)
    surface = eq.extract_flux_surface(target, npoints=n_surface_points)
    contour = eq.compute_contour_integrals(surface)
    if name == "direct":
        local = flux_surface_mercier_luc(
            surface, eq.fields, fourier_modes=fourier_modes,
            projection_iterations=projection_iterations, enforce_q=enforce_direct_q,
        )
        parameters = SimpleNamespace(
            FourierModes=fourier_modes, ProjectionIterations=projection_iterations
        )
        boundary = Boundary(
            local.theta_geo, local.R, local.Z, local.Rt, local.Zt,
            local.Rtt, local.Ztt,
        )
    else:
        options = dict(ntheta=n_fit_points, dpsi=width, n_radial_fit=n_radial_fit)
        extras = {} if fit_options is None else dict(fit_options)
        overlap = set(options) & set(extras)
        if overlap:
            raise ValueError(f"Specify {sorted(overlap)} through the preparation options.")
        options.update(extras)
        fitter, builder = fitters[name]
        parameters, boundary = fitter(eq, target, **options)
        local = builder(surface, parameters, boundary)

    source = str(volume_source).lower()
    if source == "auto":
        source = "local" if name in ("circular", "cshape", "indented-miller") else "contour"
    if source not in ("contour", "local"):
        raise ValueError("volume_source must be auto, contour, or local.")
    local_volume = local_volume_integrals(local, parameters) if source == "local" else None
    integrals = contour
    if local_volume is not None:
        V, Vprime = local_volume.V, local_volume.Vprime
        area = abs(float(trapezoid(local.R * local.Zt, local.theta_geo)))
        s = 2 * V * surface.qprime / Vprime
        alpha = -(2 * Vprime / (4 * np.pi**2)) * np.sqrt(
            V / (2 * np.pi**2 * surface.Raxis)
        ) * _MU0 * surface.pprime
        integrals = replace(
            contour, V=V, Vprime=Vprime, s=s, s_hat=s / surface.q, alpha=alpha,
            area_poloidal=area, mean_R=V / (_TWO_PI * max(area, np.finfo(float).eps)),
        )
    FFprime, pprime, qprime = (
        float(fun(target)) for fun in (eq.fields.FFprime, eq.fields.pprime, eq.fields.qprime)
    )
    Gamma, thetaPsi, thetaDiag = evaluate_theta_psi(
        local, FFprime, pprime, qprime, enforce_periodicity=enforce_periodicity
    )
    normalization = SimpleNamespace(aN=float(eq.data.rmaxis), BN=abs(float(eq.data.bcentr)))
    diagnostics = SimpleNamespace(
        qIntegralRelativeError=float(contour.q_relerr),
        qLocalRelativeError=float((local.q_check - surface.q) / max(abs(surface.q), np.finfo(float).eps)),
        FFprimePeriodicityRelativeDifference=float(
            abs(thetaDiag.FFprime_periodic - FFprime) / max(abs(FFprime), 1.0)
        ),
        periodicityClosure=float(thetaDiag.closure_used),
    )
    if hasattr(local, "diagnostics"):
        diagnostics.directGeometry = local.diagnostics
    if hasattr(parameters, "fit"):
        diagnostics.fitConverged = bool(parameters.fit.converged)
        diagnostics.fitRelativeRms = float(parameters.fit.relative_rms_distance)
    else:
        diagnostics.fitConverged = True
    geometry = SimpleNamespace(
        model=name, volumeSource=source, V=float(integrals.V),
        Vprime=float(integrals.Vprime),
        surfaceBoxKappa=float(np.ptp(local.Z) / np.ptp(local.R)),
        localVolume=local_volume,
    )
    state = SAlphaState(
        model=name, targetPsiN=target, surface=surface, integrals=integrals,
        contourIntegrals=contour, parameters=parameters, boundary=boundary,
        localEquilibrium=local, FFprime=FFprime, pprime=pprime, qprime=qprime,
        shat=float(integrals.s_hat), alpha=float(integrals.alpha),
        Gamma=Gamma, thetaPsi=thetaPsi, thetaDiag=thetaDiag,
        normalization=normalization, diagnostics=diagnostics, geometry=geometry,
        eq=eq.data, eqfunc=eq.fields, adapterDiagnostics=None,
        settings=SimpleNamespace(
            nSurfacePoints=n_surface_points, nFitPoints=n_fit_points,
            radialFitHalfWidth=width, nRadialFit=n_radial_fit,
            enforcePeriodicity=bool(enforce_periodicity), volumeSource=source,
        ),
    )
    _surface_constants(state)
    return state


def prepare_fem_salpha_equilibrium(
    equilibrium, target_psiN, *, adapter=None, adapter_nr=257, adapter_nz=None,
    **prepare_options,
):
    """Accept an adapted Equilibrium, (data, fields[, diagnostics]), or an adapter.

    A raw FEM result requires adapter(result, nr=..., nz=...) returning
    (data, fields[, diagnostics]). No unported GS solver is invoked.
    """
    _integer(adapter_nr, "adapter_nr", 4)
    if adapter_nz is not None:
        _integer(adapter_nz, "adapter_nz", 4)
    adapted = (
        adapter(equilibrium, nr=adapter_nr, nz=adapter_nz)
        if adapter is not None else equilibrium
    )
    diagnostics = None
    if hasattr(adapted, "data") and hasattr(adapted, "fields"):
        eq = adapted
    elif isinstance(adapted, (tuple, list)) and len(adapted) in (2, 3):
        eq = Equilibrium(adapted[0], adapted[1])
        if len(adapted) == 3:
            diagnostics = adapted[2]
    else:
        raise TypeError("Raw FEM results require an adapter returning (data, fields[, diagnostics]).")
    state = prepare_salpha_equilibrium(eq, target_psiN, **prepare_options)
    state.adapterDiagnostics = diagnostics
    return state


def _solve_phase(state, diagnostic, phase, options):
    def coefficients(boundary):
        return evaluate_gcf(
            diagnostic, float(phase), theta_bnd=boundary,
            a_N=state.normalization.aN, B_N=state.normalization.BN,
        )
    coeff = coefficients(options.thetaBoundary)
    # Rounded extension lengths can fall short of an arbitrary requested boundary.
    if (coeff.theta[0] > options.thetaGrid[0]
            or coeff.theta[-1] < options.thetaGrid[-1]):
        step = _TWO_PI / (len(diagnostic.theta) - 1)
        coeff = coefficients(options.thetaBoundary + step)
    if not coeff.validation.all_finite:
        raise ValueError("The ballooning coefficients are nonfinite.")
    if options.requireVectorConsistency:
        consistent = getattr(coeff.validation, "vector_consistent", None)
        if consistent is None:
            errors = [getattr(coeff.validation, name) for name in (
                "K_relative_error", "D_relative_error", "g_relative_error",
                "f_relative_error", "c_relative_error",
            )]
            consistent = np.all(np.isfinite(errors)) and max(errors) <= 1e-8
        if not consistent:
            raise ValueError(
                "Scalar/vector coefficient formulas disagree. "
                "Check the local normal and flux-sign convention."
            )
    arrays = PchipInterpolator(
        coeff.theta, np.column_stack((coeff.g, coeff.c, coeff.f)),
        axis=0, extrapolate=False,
    )(options.thetaGrid)
    if not np.all(np.isfinite(arrays)):
        raise ValueError("Coefficient interpolation did not cover the eigenvalue grid.")
    g, c, f = arrays.T
    if np.any(g <= 0) or np.any(f <= 0):
        raise ValueError("The Sturm-Liouville coefficients g and f must be positive.")
    solver = Sturm_Liouville(
        options.thetaGrid, g, c, f, bcType="dirichlet", numEigenvalues=1,
    )
    value = float(solver.eigval[0])
    if not np.isfinite(value):
        raise ValueError("The largest eigenvalue is nonfinite.")
    return value


def _solve_point(state, shat, alpha, options, *, allow_failures=False):
    shat, alpha = float(shat), float(alpha)
    if not np.all(np.isfinite([shat, alpha])):
        raise ValueError("shat and alpha must be finite.")
    qprime = float(shat_to_qprime(shat, state))
    pprime = float(alpha_to_pprime(alpha, state))
    values = np.full(options.theta0.size, np.nan)
    messages = np.full(options.theta0.size, "", dtype=object)
    used = closure = np.nan
    diagnostic = None
    try:
        _, _, diagnostic = evaluate_theta_psi(
            state.localEquilibrium, state.FFprime, pprime, qprime,
            enforce_periodicity=options.enforcePeriodicity,
        )
        used, closure = diagnostic.FFprime_used, diagnostic.closure_used
    except Exception as exc:
        if not allow_failures:
            raise
        messages[:] = f"{type(exc).__name__}: {exc}"
    if diagnostic is not None:
        for i, phase in enumerate(options.theta0):
            try:
                values[i] = _solve_phase(state, diagnostic, phase, options)
            except Exception as exc:
                if not allow_failures:
                    raise
                messages[i] = f"{type(exc).__name__}: {exc}"
    valid = np.flatnonzero(np.isfinite(values))
    best = int(valid[np.argmax(values[valid])]) if valid.size else None
    return SAlphaPoint(
        shat=shat, alpha=alpha, theta0=options.theta0.copy(),
        thetaGrid=options.thetaGrid.copy(), lambdaTheta0=values,
        lambdaMax=float(values[best]) if best is not None else np.nan,
        bestTheta0=float(options.theta0[best]) if best is not None else np.nan,
        qprime=qprime, pprime=pprime, FFprimeUsed=float(used),
        periodicityClosure=float(closure), failMessage=messages,
        validPhaseCount=int(valid.size), phaseComplete=bool(valid.size == values.size),
    )


def evaluate_salpha_point(
    state, shat, alpha, *, theta0=None, theta_boundary=5 * np.pi,
    n_balloon_points=151, enforce_periodicity=True,
    require_vector_consistency=True,
):
    """Solve all phases; raise if any phase fails, as in the MATLAB point API."""
    options = _point_options(
        theta0, theta_boundary, n_balloon_points,
        enforce_periodicity, require_vector_consistency,
    )
    _surface_constants(state)
    return _solve_point(state, shat, alpha, options)

# ----------------------------------------------------------------------
# Shared results and atomic checkpoints
# ----------------------------------------------------------------------

def _new_result(state, shat, alpha, options, *, refined=False):
    nS, nA, nT = shat.size, alpha.size, options.theta0.size
    V, Vprime, q, _ = _surface_constants(state)
    cls = RefinedSAlphaResult if refined else SAlphaScanResult
    return cls(
        model=getattr(state.localEquilibrium, "model", state.model),
        targetPsiN=state.targetPsiN, shat=shat.copy(), alpha=alpha.copy(),
        theta0=options.theta0.copy(), thetaGrid=options.thetaGrid.copy(),
        lambdaMax=np.full((nS, nA), np.nan),
        lambdaTheta0=np.full((nS, nA, nT), np.nan),
        bestTheta0=np.full((nS, nA), np.nan),
        failMessage=np.full((nS, nA, nT), "", dtype=object),
        FFprimeUsed=np.full((nS, nA), np.nan),
        periodicityClosure=np.full((nS, nA), np.nan),
        validPhaseCount=np.zeros((nS, nA), dtype=int),
        qprime=shat_to_qprime(shat, state), pprime=alpha_to_pprime(alpha, state),
        FFprimeInput=float(state.FFprime),
        equilibriumShat=float(state.shat), equilibriumAlpha=float(state.alpha),
        volume=V, Vprime=Vprime, q=q,
        geometry=getattr(state, "geometry", None),
        done=np.zeros(nS, dtype=bool), rowSeconds=np.zeros(nS),
        settings=SimpleNamespace(
            thetaBoundary=options.thetaBoundary,
            nBalloonPoints=options.nBalloonPoints,
            periodicityEnforced=options.enforcePeriodicity,
            requireVectorConsistency=options.requireVectorConsistency,
        ),
    )


def _record_point(result, row, col, point):
    result.lambdaMax[row, col] = point.lambdaMax
    result.lambdaTheta0[row, col] = point.lambdaTheta0
    result.bestTheta0[row, col] = point.bestTheta0
    result.failMessage[row, col] = point.failMessage
    result.FFprimeUsed[row, col] = point.FFprimeUsed
    result.periodicityClosure[row, col] = point.periodicityClosure
    result.validPhaseCount[row, col] = point.validPhaseCount


def _summarize(result):
    finite = result.lambdaMax[np.isfinite(result.lambdaMax)]
    result.lambdaRange = (
        np.array([finite.min(), finite.max()]) if finite.size else np.array([np.nan, np.nan])
    )
    result.hasMarginalCrossing = bool(
        finite.size and result.lambdaRange[0] <= 0 <= result.lambdaRange[1]
    )
    result.numberFailed = int(np.count_nonzero(result.failMessage != ""))
    result.phaseComplete = result.validPhaseCount == result.theta0.size
    result.numberIncompletePoints = int(np.count_nonzero(~result.phaseComplete))
    result.scanElapsedSeconds = float(np.sum(result.rowSeconds))
    if hasattr(result, "roots"):
        result.numberRoots = sum(len(roots) for roots in result.roots)
        result.elapsedSeconds = result.scanElapsedSeconds
    return result


def _hash_value(hasher, value):
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "O":
            _hash_value(hasher, value.tolist())
        else:
            hasher.update(f"A:{value.dtype.str}:{value.shape}:".encode())
            hasher.update(np.ascontiguousarray(value).tobytes())
    elif isinstance(value, np.generic):
        _hash_value(hasher, value.item())
    elif isinstance(value, dict):
        hasher.update(b"D{")
        for key in sorted(value):
            _hash_value(hasher, str(key))
            _hash_value(hasher, value[key])
        hasher.update(b"}")
    elif isinstance(value, (list, tuple)):
        hasher.update(b"L[")
        for item in value:
            _hash_value(hasher, item)
        hasher.update(b"]")
    elif hasattr(value, "__dict__"):
        _hash_value(hasher, vars(value))
    elif isinstance(value, (str, bool, float, int)) or value is None:
        encoded = json.dumps(value, ensure_ascii=True, allow_nan=True).encode()
        hasher.update(f"S:{len(encoded)}:".encode() + encoded)
    else:
        raise TypeError(f"Unsupported checkpoint signature value: {type(value).__name__}")


def _checkpoint_signature(state, shat, alpha, options, kind, tolerance=None):
    digest = hashlib.sha256()
    payload = dict(
        version=_CHECKPOINT_VERSION, kind=kind,
        localEquilibrium=state.localEquilibrium, surface=state.surface,
        integrals=state.integrals, FFprime=state.FFprime,
        normalization=state.normalization,
        shat=shat, alpha=alpha, options=vars(options), alphaTolerance=tolerance,
    )
    _hash_value(digest, payload)
    # Reject checkpoints after changes to the implemented solver code.
    for obj in (prepare_salpha_equilibrium, evaluate_gcf, Sturm_Liouville):
        try:
            filename = inspect.getsourcefile(obj)
        except (TypeError, OSError):
            filename = None
        if filename and Path(filename).is_file():
            digest.update(Path(filename).read_bytes())
    return digest.hexdigest()


_CHECKPOINT_ARRAYS = (
    "lambdaMax", "lambdaTheta0", "bestTheta0", "failMessage",
    "FFprimeUsed", "periodicityClosure", "validPhaseCount", "done", "rowSeconds",
)


def _json_value(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if hasattr(value, "__dict__"):
        return {key: _json_value(item) for key, item in vars(value).items()}
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    return value


def _save_checkpoint(filename, result, kind):
    if filename is None or str(filename) == "":
        return
    path = Path(filename).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    meta = dict(version=_CHECKPOINT_VERSION, kind=kind, signature=result.signature)
    if hasattr(result, "roots"):
        meta["roots"] = _json_value(result.roots)
    arrays = {name: getattr(result, name) for name in _CHECKPOINT_ARRAYS}
    arrays["failMessage"] = np.asarray(result.failMessage, dtype=str)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".partial")
    created = False
    try:
        # Exclusive creation propagates write errors immediately on Windows.
        with temporary.open("xb") as stream:
            created = True
            np.savez_compressed(stream, metadata=np.array(json.dumps(meta)), **arrays)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if created and temporary.exists():
            temporary.unlink()


def _restore_checkpoint(filename, result, kind, resume):
    if filename is None or str(filename) == "" or not resume:
        return result
    path = Path(filename).expanduser()
    if not path.is_file():
        return result
    with np.load(path, allow_pickle=False) as saved:
        meta = json.loads(str(saved["metadata"].item()))
        if (
            meta.get("version") != _CHECKPOINT_VERSION
            or meta.get("kind") != kind
            or meta.get("signature") != result.signature
        ):
            raise ValueError("Checkpoint settings or geometry changed; use a new file or resume=False.")
        for name in _CHECKPOINT_ARRAYS:
            if name not in saved:
                raise ValueError(f"Checkpoint array {name} is missing.")
            target = getattr(result, name)
            if saved[name].shape != target.shape:
                raise ValueError(f"Checkpoint array {name} has an incompatible shape.")
            target[...] = saved[name]
        if hasattr(result, "roots"):
            roots = meta.get("roots")
            if roots is None or len(roots) != result.shat.size:
                raise ValueError("Checkpoint roots are missing or incompatible.")
            result.roots = []
            for row in roots:
                restored = []
                for values in row:
                    values["bracket"] = np.asarray(values["bracket"], dtype=float)
                    values["lambdaBracket"] = np.asarray(values["lambdaBracket"], dtype=float)
                    restored.append(MarginalRoot(**values))
                result.roots.append(restored)
    return result


def scan_salpha_diagram(
    state, *, shat=None, alpha=None, theta0=None,
    theta_boundary=5 * np.pi, n_balloon_points=151,
    enforce_periodicity=True, require_vector_consistency=True,
    checkpoint_file=None, resume=True, verbose=True,
):
    """Scan fixed geometry, recording failed phases without clipping eigenvalues."""
    shat = _axis(np.linspace(0, 5, 31) if shat is None else shat, "shat")
    alpha = _axis(np.linspace(0, 6, 31) if alpha is None else alpha, "alpha")
    options = _point_options(
        theta0, theta_boundary, n_balloon_points,
        enforce_periodicity, require_vector_consistency,
    )
    result = _new_result(state, shat, alpha, options)
    result.signature = _checkpoint_signature(state, shat, alpha, options, "scan")
    _restore_checkpoint(checkpoint_file, result, "scan", resume)
    if verbose:
        print(f"Scanning {shat.size} shat x {alpha.size} alpha x {options.theta0.size} phases.")
    for row, s in enumerate(shat):
        if result.done[row]:
            continue
        started = time.perf_counter()
        for col, a in enumerate(alpha):
            point = _solve_point(state, s, a, options, allow_failures=True)
            _record_point(result, row, col, point)
        result.rowSeconds[row] = time.perf_counter() - started
        result.done[row] = True
        _save_checkpoint(checkpoint_file, result, "scan")
        if verbose:
            print(f"  row {row + 1}/{shat.size}: shat={s:+.5g}, "
                  f"{result.rowSeconds[row]:.1f} s")
    _summarize(result)
    if verbose:
        print(f"Finished in {result.scanElapsedSeconds:.1f} s; "
              f"{result.numberFailed} failed phases, "
              f"{result.numberIncompletePoints} incomplete points.")
    return result


def scan_fem_salpha_diagram(
    equilibrium, target_psiN, *, adapter=None, adapter_nr=257, adapter_nz=None,
    prepare_options=None, **scan_options,
):
    config = dict(model="direct", n_surface_points=512)
    if prepare_options is not None:
        config.update(dict(prepare_options))
    state = prepare_fem_salpha_equilibrium(
        equilibrium, target_psiN, adapter=adapter,
        adapter_nr=adapter_nr, adapter_nz=adapter_nz, **config,
    )
    result = scan_salpha_diagram(state, **scan_options)
    result.adapterDiagnostics = state.adapterDiagnostics
    return result, state

# ----------------------------------------------------------------------
# Marginal roots: discovery of first/second stability and bracket tracing
# ----------------------------------------------------------------------

def _refine_interval(point, left, right, tolerance, *, max_iterations=None):
    lp, rp = point(left), point(right)
    fl, fr = lp.lambdaMax, rp.lambdaMax
    direction = int(np.sign(fr - fl))
    iterations = 0
    status = "converged"
    if fl == 0:
        right, fr = left, fl
    elif fr == 0:
        left, fl = right, fr
    elif np.sign(fl) == np.sign(fr):
        raise ValueError("The refinement interval must bracket a sign change.")
    while right - left > tolerance:
        if max_iterations is not None and iterations >= max_iterations:
            status = "iteration_limit"
            break
        middle = 0.5 * (left + right)
        if middle == left or middle == right:
            status = "resolution_limit"
            break
        fm = point(middle).lambdaMax
        iterations += 1
        if fm == 0:
            left = right = middle
            fl = fr = 0.0
            break
        if np.sign(fm) == np.sign(fl):
            left, fl = middle, fm
        else:
            right, fr = middle, fm
    estimate = (
        0.5 * (left + right) if fr == fl
        else left - fl * (right - left) / (fr - fl)
    )
    estimate = float(np.clip(estimate, left, right))
    actual = point(estimate)
    return MarginalRoot(
        alpha=estimate, bracket=np.array([left, right]),
        lambdaBracket=np.array([fl, fr]), direction=direction,
        bestTheta0=actual.bestTheta0, lambdaAtRoot=actual.lambdaMax,
        iterations=iterations, status=status, evaluations=0,
    )


def scan_refined_salpha_diagram(
    state, *, shat=None, alpha=None, theta0=None,
    theta_boundary=7 * np.pi, n_balloon_points=281,
    alpha_tolerance=1e-3, checkpoint_file=None, resume=True,
    enforce_periodicity=True, require_vector_consistency=True, verbose=True,
):
    """Discover sampled sign changes and refine each using fresh eigenvalue solves.

    Features narrower than the alpha discovery spacing can be missed.
    Any failed phase stops the current row; completed checkpoint rows survive.
    """
    shat = _axis(np.linspace(0, 7, 81) if shat is None else shat,
                 "shat", nonnegative=True)
    alpha = _axis(np.linspace(0, 12, 81) if alpha is None else alpha,
                  "alpha", nonnegative=True)
    if theta0 is None:
        theta0 = _TWO_PI * (np.arange(24) - 12) / 24
    elif np.asarray(theta0).size == 0:
        raise ValueError("Refined scans require at least one explicit phase.")
    options = _point_options(
        theta0, theta_boundary, n_balloon_points,
        enforce_periodicity, require_vector_consistency,
    )
    tolerance = _positive(alpha_tolerance, "alpha_tolerance")
    result = _new_result(state, shat, alpha, options, refined=True)
    result.alphaTolerance = tolerance
    result.roots = [[] for _ in shat]
    result.signature = _checkpoint_signature(state, shat, alpha, options, "refined", tolerance)
    _restore_checkpoint(checkpoint_file, result, "refined", resume)
    for row, s in enumerate(shat):
        if result.done[row]:
            continue
        started = time.perf_counter()
        cache = {}

        def point(a):
            key = float(a)
            if key not in cache:
                cache[key] = _solve_point(state, s, key, options)
            return cache[key]

        for col, a in enumerate(alpha):
            _record_point(result, row, col, point(a))
        values = result.lambdaMax[row]
        roots = []
        # Preserve zeros at the discovery endpoints as well as interior zeros.
        for col in np.flatnonzero(values == 0):
            a = float(alpha[col])
            before = values[col - 1] if col > 0 else 0.0
            after = values[col + 1] if col + 1 < alpha.size else 0.0
            direction = 0 if np.sign(before) == np.sign(after) else int(np.sign(after - before))
            p = point(a)
            roots.append(MarginalRoot(
                alpha=a, bracket=np.array([a, a]), lambdaBracket=np.zeros(2),
                direction=direction, bestTheta0=p.bestTheta0, lambdaAtRoot=0.0,
                iterations=0, status="sampled_zero", evaluations=0,
            ))
        for col in range(alpha.size - 1):
            fl, fr = values[col], values[col + 1]
            if fl == 0 or fr == 0 or np.sign(fl) == np.sign(fr):
                continue
            previous_count = len(cache)
            root = _refine_interval(
                point, float(alpha[col]), float(alpha[col + 1]), tolerance
            )
            root.evaluations = len(cache) - previous_count
            roots.append(root)
        result.roots[row] = sorted(roots, key=lambda item: item.alpha)
        result.rowSeconds[row] = time.perf_counter() - started
        result.done[row] = True
        _save_checkpoint(checkpoint_file, result, "refined")
        if verbose:
            print(f"  refined row {row + 1}/{shat.size}: shat={s:.5g}, "
                  f"{len(roots)} boundaries, {result.rowSeconds[row]:.1f} s")
    return _summarize(result)


def trace_salpha_marginal(
    state, shat, *, alpha_bracket=(0.0, 3.0), alpha_tolerance=5e-3,
    max_iterations=20, theta0=None, theta_boundary=5 * np.pi,
    n_balloon_points=151, enforce_periodicity=True,
    require_vector_consistency=True, verbose=True,
):
    """Refine one stable-left/unstable-right bracket per shat.

    Use scan_refined_salpha_diagram to discover multiple stability regions.
    lambdaAtMarginal is evaluated at marginalAlpha, not at a nearby endpoint.
    """
    shat = _vector(shat, "shat")
    bracket = _axis(alpha_bracket, "alpha_bracket")
    if bracket.size != 2:
        raise ValueError("alpha_bracket must contain exactly two values.")
    tolerance = _positive(alpha_tolerance, "alpha_tolerance")
    _integer(max_iterations, "max_iterations", 1)
    options = _point_options(
        theta0, theta_boundary, n_balloon_points,
        enforce_periodicity, require_vector_consistency,
    )
    _surface_constants(state)
    count = shat.size
    result = MarginalTrace(
        shat=shat.copy(), marginalAlpha=np.full(count, np.nan),
        lambdaAtMarginal=np.full(count, np.nan), bestTheta0=np.full(count, np.nan),
        bracketWidth=np.full(count, np.nan), iterationCount=np.zeros(count, dtype=int),
        evaluationCount=np.zeros(count, dtype=int),
        status=np.full(count, "", dtype=object), message=np.full(count, "", dtype=object),
        alphaBracket=bracket.copy(), alphaTolerance=tolerance,
    )
    for row, s in enumerate(shat):
        cache = {}

        def point(a):
            key = float(a)
            if key not in cache:
                result.evaluationCount[row] += 1
                cache[key] = _solve_point(state, s, key, options)
            return cache[key]

        try:
            left, right = map(float, bracket)
            lp, rp = point(left), point(right)
            if lp.lambdaMax > 0 or rp.lambdaMax < 0:
                result.status[row] = "unbracketed"
                result.message[row] = (
                    f"lambda(left)={lp.lambdaMax:+.6e}, lambda(right)={rp.lambdaMax:+.6e}"
                )
            else:
                root = _refine_interval(
                    point, left, right, tolerance, max_iterations=max_iterations
                )
                result.marginalAlpha[row] = root.alpha
                result.lambdaAtMarginal[row] = root.lambdaAtRoot
                result.bestTheta0[row] = root.bestTheta0
                result.bracketWidth[row] = np.ptp(root.bracket)
                result.iterationCount[row] = root.iterations
                result.status[row] = root.status
        except Exception as exc:
            result.status[row] = "failed"
            result.message[row] = f"{type(exc).__name__}: {exc}"
        if verbose:
            print(f"  marginal row {row + 1}/{count}: shat={s:+.5g}, "
                  f"alpha={result.marginalAlpha[row]:+.6f}, {result.status[row]}")
    result.numberConverged = int(np.count_nonzero(result.status == "converged"))
    result.numberUnbracketed = int(np.count_nonzero(result.status == "unbracketed"))
    result.numberFailed = int(np.count_nonzero(result.status == "failed"))
    result.numberLimited = int(np.count_nonzero(np.isin(
        result.status, ("iteration_limit", "resolution_limit")
    )))
    return result


# ----------------------------------------------------------------------
# Plotting (matplotlib is imported only when this function is called)
# ----------------------------------------------------------------------

def _has_adjacent_crossing(values):
    for left, right in ((values[:, :-1], values[:, 1:]), (values[:-1], values[1:])):
        valid = np.isfinite(left) & np.isfinite(right)
        if np.any(valid & (((left < 0) & (right >= 0)) | ((left > 0) & (right <= 0)))):
            return True
    return False


def plot_salpha_diagram(
    result, *, show=False, show_phase=True, title=None,
    contour_levels=30, mask_incomplete=True,
):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    _integer(contour_levels, "contour_levels", 2)
    shat, alpha = _axis(result.shat, "shat"), _axis(result.alpha, "alpha")
    values = np.array(result.lambdaMax, dtype=float, copy=True)
    phases = np.array(result.bestTheta0, dtype=float, copy=True)
    if values.shape != (shat.size, alpha.size) or phases.shape != values.shape:
        raise ValueError("lambdaMax and bestTheta0 must have shape (n_shat, n_alpha).")
    if mask_incomplete and hasattr(result, "phaseComplete"):
        values[~result.phaseComplete] = np.nan
        phases[~result.phaseComplete] = np.nan
    stability, ax = plt.subplots()
    stability.canvas.manager.set_window_title("s-alpha stability diagram")
    finite = np.isfinite(values)
    if np.any(finite):
        image = ax.contourf(
            alpha, shat, np.ma.masked_invalid(values),
            levels=contour_levels, cmap="turbo",
        )
        stability.colorbar(image, ax=ax, label=r"$\lambda_{\max}$")
        if _has_adjacent_crossing(values):
            ax.contour(alpha, shat, values, levels=[0.0], colors="black", linewidths=2.2)
    else:
        ax.text(0.5, 0.5, "No valid eigenvalues", ha="center", va="center", transform=ax.transAxes)
        ax.set(xlim=(alpha[0], alpha[-1]), ylim=(shat[0], shat[-1]))
    eq_a, eq_s = result.equilibriumAlpha, result.equilibriumShat
    if np.isfinite(eq_a) and np.isfinite(eq_s) and alpha[0] <= eq_a <= alpha[-1] and shat[0] <= eq_s <= shat[-1]:
        ax.plot(eq_a, eq_s, marker="*", markersize=14,
                markerfacecolor="gold", markeredgecolor="black", linestyle="none")
    label = getattr(result, "model", "local equilibrium").replace("-", " ")
    ax.set(
        xlabel=r"$\alpha$", ylabel=r"$\hat{s}$",
        title=title if title else rf"{label}: $\psi_N={result.targetPsiN:.3f}$",
    )
    ax.grid(True, alpha=0.25)
    stability.tight_layout()
    figures = SimpleNamespace(stability=stability, stabilityAxes=ax, phase=None, phaseAxes=None)
    if show_phase:
        phase_figure, phase_ax = plt.subplots()
        phase_figure.canvas.manager.set_window_title("Phase of maximum ballooning eigenvalue")
        grid_a, grid_s = np.meshgrid(alpha, shat)
        valid = np.isfinite(phases)
        sampled = np.sort(_vector(result.theta0, "theta0")) / np.pi
        if sampled.size == 1:
            edges = np.array([sampled[0] - 1, sampled[0] + 1])
        else:
            middle = 0.5 * (sampled[:-1] + sampled[1:])
            edges = np.r_[sampled[0] - (middle[0] - sampled[0]),
                          middle, sampled[-1] + (sampled[-1] - middle[-1])]
        hue = np.mod((sampled + 1) / 2, 1)
        colors = ListedColormap(plt.get_cmap("hsv")(hue))
        norm = BoundaryNorm(edges, sampled.size)
        scatter = phase_ax.scatter(
            grid_a[valid], grid_s[valid], c=phases[valid] / np.pi,
            cmap=colors, norm=norm, marker="s", s=58,
            edgecolors=(0.18, 0.18, 0.18), linewidths=0.25,
        )
        ticks = sampled[::max(1, int(np.ceil(sampled.size / 12)))]
        bar = phase_figure.colorbar(scatter, ax=phase_ax, ticks=ticks)
        bar.set_label(r"$\theta_{0,\max}/\pi$")
        phase_ax.set(
            xlabel=r"$\alpha$", ylabel=r"$\hat{s}$",
            title="Phase of maximum ballooning eigenvalue",
            xlim=(alpha[0] - 0.05 * np.ptp(alpha), alpha[-1] + 0.05 * np.ptp(alpha)),
            ylim=(shat[0] - 0.05 * np.ptp(shat), shat[-1] + 0.05 * np.ptp(shat)),
        )
        phase_ax.grid(True, alpha=0.25)
        phase_figure.tight_layout()
        figures.phase, figures.phaseAxes = phase_figure, phase_ax
    if show:
        plt.show()
    return figures


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geqdsk")
    parser.add_argument("--psiN", type=float, default=0.7)
    parser.add_argument("--model", choices=("miller", "circular", "cshape", "indented", "direct"),
                        default="miller")
    parser.add_argument("--nshat", type=int, default=11)
    parser.add_argument("--nalpha", type=int, default=11)
    parser.add_argument("--ntheta0", type=int, default=8)
    parser.add_argument("--shat-end", type=float, default=5.0)
    parser.add_argument("--alpha-end", type=float, default=6.0)
    parser.add_argument("--theta-boundary", type=float, default=5 * np.pi)
    parser.add_argument("--n-balloon-points", type=int, default=151)
    parser.add_argument("--surface-points", type=int, default=512)
    parser.add_argument("--fit-points", type=int, default=256)
    parser.add_argument("--fit-iterations", type=int, default=200)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    _integer(args.ntheta0, "ntheta0", 1)
    fit_config = (
        {"max_iterations": args.fit_iterations} if args.model in ("cshape", "indented") else None
    )
    state = prepare_salpha_equilibrium(
        args.geqdsk, args.psiN, model=args.model,
        n_surface_points=args.surface_points, n_fit_points=args.fit_points,
        fit_options=fit_config,
    )
    if not state.diagnostics.fitConverged:
        warnings.warn("The nonlinear geometry fit has not converged; inspect state.diagnostics.",
                      RuntimeWarning)
    phases = _TWO_PI * (np.arange(args.ntheta0) - args.ntheta0 // 2) / args.ntheta0
    result = scan_salpha_diagram(
        state, shat=np.linspace(0, args.shat_end, args.nshat),
        alpha=np.linspace(0, args.alpha_end, args.nalpha),
        theta0=phases, theta_boundary=args.theta_boundary,
        n_balloon_points=args.n_balloon_points, checkpoint_file=args.checkpoint,
    )
    print(f"lambda range: {result.lambdaRange}")
    if args.plot:
        plot_salpha_diagram(result, show=True)
