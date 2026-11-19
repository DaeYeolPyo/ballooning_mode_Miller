"""Normalized-flux profiles for the fixed-boundary P3 GS solver.

pprime  means dp/dpsi_N.
ffprime means F*dF/dpsi_N = 0.5*d(F**2)/dpsi_N.

No conversion by the dimensional flux span is performed here.
Pass these arrays directly to fem.assemble_rhs(), which divides by dpsi.

The reconstruction routines accept nonuniform profile grids.
The current FEM RHS requires a uniform grid from 0 to 1.
"""

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import PchipInterpolator


MU0 = 4.0 * np.pi * 1.0e-7
_REPRESENTATION = "direct-FFprime-reduced-basis"

__all__ = [
    "FProfile", "FFprimeModel", "FFprimeOptions", "FFprimeInitial",
    "InitializationInfo", "CylindricalInfo", "QTarget",
    "validate_profile_grid", "construct_ffprime_basis",
    "reconstruct_g_from_ffprime", "reconstruct_pressure",
    "evaluate_ffprime_model", "initialize_ffprime_model",
    "evaluate_q_target", "normalize_q_target",
]


def _vector(values, name, size=None, allow_scalar=False):
    a = np.asarray(values)
    if np.iscomplexobj(a):
        raise ValueError(f"{name} must be real.")
    a = np.array(a, dtype=float, copy=True)
    if allow_scalar and a.ndim == 0:
        a = a.reshape(1)
    if a.ndim not in (1, 2) or (a.ndim == 2 and 1 not in a.shape):
        raise ValueError(f"{name} must be a vector.")
    a = a.reshape(-1)
    if size is not None and a.size != size:
        raise ValueError(f"{name} must contain {size} values.")
    if not a.size or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must contain finite values.")
    return a


def _scalar(value, name, nonzero=False, positive=False):
    a = np.asarray(value)
    if a.ndim != 0 or np.iscomplexobj(a) or a.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a real scalar.")
    value = float(a)
    if (
        not np.isfinite(value)
        or (nonzero and value == 0)
        or (positive and value <= 0)
    ):
        raise ValueError(f"{name} has an invalid value.")
    return value


def _integer(value, name, lower=1, upper=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer.")
    value = int(value)
    if value < lower or (upper is not None and value > upper):
        raise ValueError(f"{name} is outside the supported range.")
    return value


def _boolean(value, name):
    if not isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be boolean.")
    return bool(value)


def _readonly(a):
    a.setflags(write=False)
    return a


def _sample_grid(values, name):
    s = _vector(values, name)
    if s.size < 2 or np.any(np.diff(s) <= 0):
        raise ValueError(f"{name} must have at least two strictly increasing samples.")
    return s


def validate_profile_grid(psi_n, *, require_uniform=False):
    """Return a copied, canonical grid spanning psi_N=0 to 1."""
    s = _sample_grid(psi_n, "psi_n")
    tolerance = 100.0 * np.spacing(max(1.0, float(np.max(np.abs(s)))))
    if abs(s[0]) > tolerance or abs(s[-1] - 1.0) > tolerance:
        raise ValueError("psi_n must span the magnetic axis 0 and boundary 1.")
    s[0], s[-1] = 0.0, 1.0
    if np.any(s < 0) or np.any(s > 1) or np.any(np.diff(s) <= 0):
        raise ValueError("Invalid normalized-flux grid after endpoint normalization.")
    if require_uniform:
        uniform = np.linspace(0.0, 1.0, s.size)
        if np.max(np.abs(s - uniform)) > 100.0 * np.spacing(float(s.size)):
            raise ValueError("The FEM RHS requires a uniform psi_n grid.")
    return s


@dataclass(frozen=True)
class FProfile:
    psi_n: np.ndarray
    g: np.ndarray
    f: np.ndarray
    ffprime: np.ndarray
    dg_dpsi_n: np.ndarray
    f_boundary: float
    integrated_2ffprime: np.ndarray
    basis: Optional[np.ndarray] = field(default=None, repr=False)
    coefficients: Optional[np.ndarray] = None
    representation: str = "sampled-FFprime"

    @property
    def g_boundary(self):
        return self.f_boundary**2

    @property
    def sign_f(self):
        return int(np.sign(self.f_boundary))

    @property
    def minimum_g(self):
        return float(np.min(self.g))

    @property
    def minimum_g_index(self):
        """Zero-based sample index."""
        return int(np.argmin(self.g))

    @property
    def boundary_error(self):
        return float(abs(self.g[-1] - self.g_boundary))


@dataclass(frozen=True)
class CylindricalInfo:
    leading_coefficient: float
    coefficient_history: np.ndarray
    flux_span_magnitude_history: np.ndarray
    final_cylinder_ffprime: np.ndarray


@dataclass(frozen=True)
class InitializationInfo:
    kind: str
    source_ffprime: np.ndarray
    cylindrical: Optional[CylindricalInfo] = None


@dataclass(frozen=True)
class FFprimeModel:
    psi_n: np.ndarray
    coefficients: np.ndarray
    f_boundary: float
    enforce_zero_edge: bool = True
    initialization: Optional[InitializationInfo] = field(default=None, repr=False)
    normalization: str = field(default="fixed_fb", init=False)
    representation: str = field(default=_REPRESENTATION, init=False)

    def __post_init__(self):
        s = validate_profile_grid(self.psi_n)
        c = _vector(self.coefficients, "coefficients")
        _integer(c.size, "n_basis", upper=12)
        fb = _scalar(self.f_boundary, "f_boundary", nonzero=True)
        edge = _boolean(self.enforce_zero_edge, "enforce_zero_edge")
        object.__setattr__(self, "psi_n", _readonly(s))
        object.__setattr__(self, "coefficients", _readonly(c))
        object.__setattr__(self, "f_boundary", fb)
        object.__setattr__(self, "enforce_zero_edge", edge)

    @property
    def n_basis(self):
        return self.coefficients.size

    @property
    def sign_f(self):
        return int(np.sign(self.f_boundary))

    def with_coefficients(self, coefficients):
        """Return an independent model for an outer-iteration trial."""
        return replace(self, coefficients=_vector(coefficients, "coefficients", self.n_basis))


@dataclass(frozen=True)
class FFprimeOptions:
    n_basis: int = 3
    enforce_zero_edge: bool = True
    cylinder_iterations: int = 6
    cylinder_relaxation: float = 0.5
    cylinder_maximum_g_drop: float = 0.8

    def __post_init__(self):
        object.__setattr__(self, "n_basis", _integer(self.n_basis, "n_basis", upper=12))
        object.__setattr__(self, "enforce_zero_edge", _boolean(self.enforce_zero_edge, "enforce_zero_edge"))
        object.__setattr__(self, "cylinder_iterations", _integer(self.cylinder_iterations, "cylinder_iterations"))
        relaxation = _scalar(self.cylinder_relaxation, "cylinder_relaxation")
        maximum = _scalar(self.cylinder_maximum_g_drop, "cylinder_maximum_g_drop")
        if not 0 < relaxation <= 1 or not 0 < maximum < 1:
            raise ValueError("Require 0 < relaxation <= 1 and 0 < maximum_g_drop < 1.")
        object.__setattr__(self, "cylinder_relaxation", relaxation)
        object.__setattr__(self, "cylinder_maximum_g_drop", maximum)


@dataclass(frozen=True)
class FFprimeInitial:
    kind: str = "cylindrical"
    psi_n: Optional[np.ndarray] = None
    ffprime: Optional[np.ndarray] = None

    def __post_init__(self):
        if not isinstance(self.kind, str):
            raise ValueError("Initialization kind must be text.")
        kind = self.kind.strip().lower()
        kind = {"zero": "constant", "constant_g": "constant"}.get(kind, kind)
        if kind not in ("constant", "ffprime", "cylindrical"):
            raise ValueError("Initialization kind must be constant, ffprime or cylindrical.")
        object.__setattr__(self, "kind", kind)
        if self.ffprime is not None:
            v = _vector(self.ffprime, "initial.ffprime")
            if v.size < 2:
                raise ValueError("initial.ffprime needs at least two samples.")
            object.__setattr__(self, "ffprime", _readonly(v))
        if self.psi_n is not None:
            s = _sample_grid(self.psi_n, "initial.psi_n")
            if self.ffprime is not None and s.size != self.ffprime.size:
                raise ValueError("Initial grid and ffprime must have equal lengths.")
            object.__setattr__(self, "psi_n", _readonly(s))
        if kind == "ffprime" and self.ffprime is None:
            raise ValueError("The ffprime initialization requires initial.ffprime.")


@dataclass(frozen=True)
class QTarget:
    """Sampled prescribed q with MATLAB-compatible PCHIP extrapolation."""
    psi_n: np.ndarray
    q: np.ndarray
    _interpolator: object = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        s = _sample_grid(self.psi_n, "q_target.psi_n")
        q = _vector(self.q, "q_target.q", s.size)
        if np.any(q == 0):
            raise ValueError("q target samples must be nonzero.")
        object.__setattr__(self, "psi_n", _readonly(s))
        object.__setattr__(self, "q", _readonly(q))
        object.__setattr__(
            self, "_interpolator", PchipInterpolator(s, q, extrapolate=True)
        )

    def __call__(self, psi_n):
        s = _vector(psi_n, "query psi_n", allow_scalar=True)
        values = np.asarray(self._interpolator(s))
        if not np.all(np.isfinite(values)) or np.any(values == 0):
            raise ValueError("Evaluated q must be finite and nonzero.")
        return values


def construct_ffprime_basis(psi_n, n_basis, enforce_zero_edge=True):
    """Return basis[:, j] = (1-s)*P_j(2*s-1), j=0..n_basis-1.

    With enforce_zero_edge=False, return just P_j(2*s-1).
    Query points need not be ordered or include the endpoints.
    """
    s = _vector(psi_n, "psi_n", allow_scalar=True)
    n = _integer(n_basis, "n_basis", upper=12)
    edge = _boolean(enforce_zero_edge, "enforce_zero_edge")
    if np.any(s < 0) or np.any(s > 1):
        raise ValueError("Basis query points must lie in [0, 1].")
    x = 2.0 * s - 1.0
    basis = np.empty((s.size, n))
    basis[:, 0] = 1.0
    if n > 1:
        basis[:, 1] = x
    for degree in range(2, n):
        basis[:, degree] = (
            (2 * degree - 1) * x * basis[:, degree - 1]
            - (degree - 1) * basis[:, degree - 2]
        ) / degree
    if edge:
        basis *= (1.0 - s)[:, None]
    return basis


def reconstruct_g_from_ffprime(psi_n, ffprime, f_boundary):
    """Reconstruct G=F**2 using the original cumulative trapezoidal rule.

    G(s) = Fb**2 + integral_0^s 2*ffprime - integral_0^1 2*ffprime.
    F keeps the constant sign of Fb; nonpositive G is rejected.
    """
    s = validate_profile_grid(psi_n)
    ff = _vector(ffprime, "ffprime", s.size)
    fb = _scalar(f_boundary, "f_boundary", nonzero=True)
    derivative = 2.0 * ff
    integrated = cumulative_trapezoid(derivative, s, initial=0.0)
    gb = fb * fb
    if not np.isfinite(gb):
        raise ValueError("f_boundary squared must be finite.")
    g = gb + integrated - integrated[-1]
    if not np.all(np.isfinite(g)):
        raise ValueError("The reconstructed F^2 profile contains NaN or Inf.")
    minimum_index = int(np.argmin(g))
    if g[minimum_index] <= 0:
        raise ValueError(
            f"Nonpositive F^2: min(G)={g[minimum_index]:.6e} "
            f"at psi_n={s[minimum_index]:.6g}."
        )
    f = np.sign(fb) * np.sqrt(g)
    scale = max(gb, float(np.max(np.abs(g))), 1.0)
    if abs(g[-1] - gb) > 100.0 * np.spacing(scale):
        raise ValueError("Reconstructed G does not preserve the boundary value.")
    g[-1], f[-1] = gb, fb
    return FProfile(
        _readonly(s), _readonly(g), _readonly(f), _readonly(ff),
        _readonly(derivative), fb, _readonly(integrated),
    )


def reconstruct_pressure(psi_n, pprime, p_edge=0.0):
    """Return p(s) = p_edge + integral_1^s (dp/dpsi_N) ds."""
    s = validate_profile_grid(psi_n)
    pp = _vector(pprime, "pprime", s.size)
    edge = _scalar(p_edge, "p_edge")
    integrated = cumulative_trapezoid(pp, s, initial=0.0)
    pressure = edge + integrated - integrated[-1]
    if not np.all(np.isfinite(pressure)):
        raise ValueError("The reconstructed pressure contains NaN or Inf.")
    pressure[-1] = edge
    return _readonly(pressure)


def evaluate_ffprime_model(model, psi_n=None):
    """Evaluate coefficients, then reconstruct the signed F profile."""
    if not isinstance(model, FFprimeModel):
        raise TypeError("model must be an FFprimeModel.")
    s = validate_profile_grid(model.psi_n if psi_n is None else psi_n)
    basis = construct_ffprime_basis(s, model.n_basis, model.enforce_zero_edge)
    ff = basis @ model.coefficients
    profile = reconstruct_g_from_ffprime(s, ff, model.f_boundary)
    return replace(
        profile, basis=_readonly(basis),
        coefficients=_readonly(model.coefficients.copy()),
        representation=_REPRESENTATION,
    )


def evaluate_q_target(q_target, psi_n):
    """Evaluate a scalar, callable, sampled target or same-length vector."""
    s = _vector(psi_n, "query psi_n", allow_scalar=True)
    if isinstance(q_target, Mapping):
        q_target = QTarget(**q_target)
    if isinstance(q_target, QTarget) or callable(q_target):
        values = q_target(s)
    else:
        a = np.asarray(q_target)
        if a.ndim == 0:
            values = np.full(s.size, _scalar(q_target, "q_target", nonzero=True))
        else:
            values = q_target
    values = _vector(values, "evaluated q_target", s.size)
    if np.any(values == 0):
        raise ValueError("Evaluated q must be nonzero.")
    return values


def normalize_q_target(q_target, psi_n):
    """Bind a numeric q vector to the common profile grid.

    Call before evaluating on the different contour grid used by mode 3.
    """
    s = validate_profile_grid(psi_n)
    if isinstance(q_target, Mapping):
        return QTarget(**q_target)
    if isinstance(q_target, QTarget) or callable(q_target):
        return q_target
    a = np.asarray(q_target)
    if a.ndim == 0:
        return _scalar(q_target, "q_target", nonzero=True)
    return QTarget(s, _vector(q_target, "q_target", s.size))


def _coerce_options(options):
    if options is None:
        return FFprimeOptions()
    if isinstance(options, Mapping):
        return FFprimeOptions(**options)
    if not isinstance(options, FFprimeOptions):
        raise TypeError("options must be FFprimeOptions or a mapping.")
    return options


def _coerce_initial(initial):
    if initial is None:
        return FFprimeInitial()
    if isinstance(initial, str):
        return FFprimeInitial(kind=initial)
    if isinstance(initial, Mapping):
        return FFprimeInitial(**initial)
    if not isinstance(initial, FFprimeInitial):
        raise TypeError("initial must be FFprimeInitial, text or a mapping.")
    return initial


def _estimate_cylindrical(s, pprime, q, fb, r0, a_minor, options):
    """Cylindrical leading-(1-s) estimate with one-dimensional sources.

    The MATLAB source transposes ppval at a column grid, creating an
    unintended N-by-N source. This implementation preserves the intended
    elementwise formula with one derivative and source value per sample.
    """
    q_magnitude = np.abs(q)
    fb_magnitude = abs(fb)
    f = np.full(s.size, fb_magnitude)
    coefficient = 0.0
    leading_basis = 1.0 - s
    fit_mask = s <= 0.9
    history = np.empty(options.cylinder_iterations)
    span_history = np.empty_like(history)

    for iteration in range(options.cylinder_iterations):
        q_over_f_integral = trapezoid(q_magnitude / f, s)
        if not np.isfinite(q_over_f_integral) or q_over_f_integral <= 0:
            raise ValueError("Failed to estimate the cylindrical flux span.")
        span = a_minor**2 / (2.0 * r0 * q_over_f_integral)
        radius_squared = 2.0 * r0 * span * cumulative_trapezoid(
            q_magnitude / f, s, initial=0.0
        )
        if not np.all(np.isfinite(radius_squared)) or radius_squared[-1] <= 0:
            raise ValueError("Failed to reconstruct the cylindrical radius.")
        radius_squared *= a_minor**2 / radius_squared[-1]
        flux_function = radius_squared * f / (r0**2 * q_magnitude)
        derivative = PchipInterpolator(s, flux_function).derivative()(s)
        cylinder_ffprime = (
            -f / q_magnitude * derivative - MU0 * r0**2 * pprime
        )
        numerator = leading_basis[fit_mask] @ cylinder_ffprime[fit_mask]
        denominator = leading_basis[fit_mask] @ leading_basis[fit_mask]
        candidate = numerator / max(denominator, np.finfo(float).eps)
        candidate = min(candidate, options.cylinder_maximum_g_drop * fb_magnitude**2)
        candidate = max(candidate, -fb_magnitude**2)
        omega = options.cylinder_relaxation
        coefficient = (1.0 - omega) * coefficient + omega * candidate
        reconstructed = reconstruct_g_from_ffprime(s, coefficient * leading_basis, fb)
        f = np.abs(reconstructed.f)
        history[iteration] = coefficient
        span_history[iteration] = span

    if not np.all(np.isfinite(cylinder_ffprime)):
        raise ValueError("The cylindrical source contains NaN or Inf.")
    info = CylindricalInfo(
        coefficient, _readonly(history), _readonly(span_history),
        _readonly(cylinder_ffprime.copy()),
    )
    return coefficient, cylinder_ffprime, info


def initialize_ffprime_model(
    psi_n, pprime, q_target=None, *, r0, b0, a_minor,
    options=None, initial=None, f_boundary=None,
    normalization="fixed_fb", f_sign=None,
):
    """Initialize the mode-3 direct-FFprime predictor.

    This interface takes grid and dimensional data explicitly, avoiding
    a dependency on the later GS problem/solver classes.

    options:
        FFprimeOptions or a mapping with the same snake_case keys.
    initial:
        "constant", "cylindrical", or FFprimeInitial(kind="ffprime", ...).
        A mapping with FFprimeInitial's field names is also accepted.
    f_boundary:
        Defaults to r0*b0. Only fixed_fb normalization is supported.

    The cylindrical estimate supplies an initial guess only; subsequent
    mode-3 updates must be checked with a complete two-dimensional GS solve.
    """
    s = validate_profile_grid(psi_n)
    pp = _vector(pprime, "pprime", s.size)
    r0 = _scalar(r0, "r0", positive=True)
    b0 = _scalar(b0, "b0")
    a_minor = _scalar(a_minor, "a_minor", positive=True)
    options = _coerce_options(options)
    initial = _coerce_initial(initial)

    if not isinstance(normalization, str):
        raise ValueError("normalization must be fixed_fb.")
    mode = normalization.strip().lower().replace("-", "_")
    if mode == "fixedfb":
        mode = "fixed_fb"
    if mode != "fixed_fb":
        raise ValueError("Only fixed_fb normalization is supported.")
    fb = _scalar(r0 * b0 if f_boundary is None else f_boundary, "f_boundary", nonzero=True)
    if f_sign is not None:
        sign = _scalar(f_sign, "f_sign")
        if sign not in (-1.0, 1.0) or sign != np.sign(fb):
            raise ValueError("f_sign must be +/-1 and agree with f_boundary.")

    basis = construct_ffprime_basis(s, options.n_basis, options.enforce_zero_edge)
    coefficients = np.zeros(options.n_basis)
    source = np.zeros_like(s)
    cylinder_info = None

    if initial.kind == "ffprime":
        values = initial.ffprime
        if initial.psi_n is None:
            if values.size != s.size:
                raise ValueError("initial.psi_n is required for a different sample count.")
            initial_grid = s
        else:
            initial_grid = initial.psi_n
        # interp1(..., 'pchip') extrapolates by default in MATLAB.
        source = PchipInterpolator(initial_grid, values, extrapolate=True)(s)
        normal = basis.T @ basis + 1.0e-12 * np.eye(options.n_basis)
        coefficients = np.linalg.solve(normal, basis.T @ source)
    elif initial.kind == "cylindrical":
        q = evaluate_q_target(q_target, s)
        coefficients[0], source, cylinder_info = _estimate_cylindrical(
            s, pp, q, fb, r0, a_minor, options
        )

    model = FFprimeModel(
        s, coefficients, fb, options.enforce_zero_edge,
        InitializationInfo(initial.kind, _readonly(source), cylinder_info),
    )
    # Check physical admissibility of the actual basis representation.
    evaluate_ffprime_model(model)
    return model
