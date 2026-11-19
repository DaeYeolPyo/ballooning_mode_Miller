"""Axisymmetric DESC H5 geometry and controlled interpolation experiments.

The source equilibrium is never changed. Fourier-Zernike conventions follow
DESC v0.17.3 basis.py; Jacobi evaluation uses SciPy, without JAX or DESC imports.
Only axisymmetric FourierZernikeBasis / HermiteSplineProfile H5 files are accepted.
"""
from pathlib import Path
import hashlib
import json

import h5py
import numpy as np
from scipy.constants import mu_0
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import CubicHermiteSpline, CubicSpline, PPoly, PchipInterpolator
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import eval_jacobi

TWO_PI = 2 * np.pi
DERIVATIVES = {"": (0, 0), "_r": (1, 0), "_t": (0, 1),
               "_rr": (2, 0), "_rt": (1, 1), "_tt": (0, 2)}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def radial(l, m, r, dr=0):
    """R_l^m = r^m P_((l-m)/2)^(0,m)(2r^2-1), including axis limits."""
    m = abs(int(m))
    n = (int(l)-m)//2
    if n < 0 or (int(l)-m) % 2:
        raise ValueError("Invalid Zernike mode")
    x = 2*np.asarray(r)**2-1
    p = eval_jacobi(n, 0, m, x)
    if dr == 0:
        return r**m*p
    dp = .5*(n+m+1)*eval_jacobi(n-1, 1, m+1, x) if n else np.zeros_like(x)
    if dr == 1:
        return (m*r**max(m-1, 0)*p if m else 0) + 4*r**(m+1)*dp
    ddp = .25*(n+m+1)*(n+m+2)*eval_jacobi(n-2, 2, m+2, x) if n > 1 else np.zeros_like(x)
    return ((m*(m-1)*r**max(m-2, 0)*p if m >= 2 else 0)
            + 4*(2*m+1)*r**m*dp + 16*r**(m+2)*ddp)


def trig(m, theta, dt=0):
    m = np.asarray(m)
    return abs(m)**dt * np.sin(np.asarray(theta)[:, None]*abs(m)
                              + (m >= 0)*np.pi/2 + dt*np.pi/2)


def invert_angle(evaluate, target):
    """Solve t+lambda(r,t)=theta_PEST on the unwrapped real line."""
    t = np.asarray(target, dtype=float).copy()
    for _ in range(70):
        value, slope = evaluate(t)
        if np.any(1+slope <= 0):
            raise ValueError("Nonmonotone straight-field-line angle")
        residual = t+value-target
        if np.max(abs(residual)) < 3e-12:
            return t
        t -= np.clip(residual/(1+slope), -.5, .5)
    raise ValueError("Straight-field-line angle inversion did not converge")


def to_pest(raw):
    """Analytic coordinate chain rule at fixed DESC toroidal-flux radius rho."""
    d = 1+raw["lambda_t"]
    ar = -raw["lambda_r"]/d
    arr = -(raw["lambda_rr"]+2*raw["lambda_rt"]*ar+raw["lambda_tt"]*ar**2)/d
    at = 1/d
    att = -raw["lambda_tt"]/d**3
    art = -(raw["lambda_rt"]+raw["lambda_tt"]*ar)/d**2
    out = {}
    for k in ("R", "Z"):
        v, vr, vt, vrr, vrt, vtt = (raw[k+s] for s in DERIVATIVES)
        out.update({k: v, k+"_r": vr+vt*ar, k+"_t": vt*at,
                    k+"_rr": vrr+2*vrt*ar+vtt*ar**2+vt*arr,
                    k+"_rt": (vrt+vtt*ar)*at+vt*art,
                    k+"_tt": vtt*at**2+vt*att})
    for suffix in DERIVATIVES:
        out["lambda"+suffix] = np.zeros_like(raw["R"])
    return out


class PeriodicFourier:
    """Trigonometric interpolant of uniformly spaced, real periodic samples."""
    def __init__(self, values):
        self.n = len(values)-1
        self.k = np.arange(self.n//2+1)
        self.coeff = np.fft.rfft(values[:-1])/self.n
        self.weight = np.full(len(self.k), 2.)
        self.weight[0] = 1.
        if self.n % 2 == 0:
            self.weight[-1] = 1.  # a real Nyquist cosine has no distinct conjugate

    def __call__(self, theta, nu=0):
        return np.real(np.exp(1j*np.asarray(theta)[:, None]*self.k)
                       @ (self.weight*self.coeff*(1j*self.k)**nu))


class H5Equilibrium:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.sha = digest(path)
        self.modes, self.coeffs = {}, {}
        with h5py.File(path, "r") as h:
            if int(h["_N"][()]) != 0:
                raise ValueError("This diagnostic requires an axisymmetric equilibrium")
            self.resolution = [int(h[k][()]) for k in ("_L", "_M", "_N")]
            self.version = h["__version__"][()].decode()
            self.Psi = float(h["_Psi"][()])
            if self.Psi <= 0:
                raise ValueError("Require positive total toroidal flux")
            for key, store in (("R", "R"), ("Z", "Z"), ("lambda", "L")):
                if h[f"_{store}_basis/__class__"][()].decode().split(".")[-1] != "FourierZernikeBasis":
                    raise ValueError("Unsupported spectral basis")
                modes = h[f"_{store}_basis/_modes"][()]
                if np.any(modes[:, 2]):
                    raise ValueError("Toroidal Fourier modes are unsupported")
                self.modes[key], self.coeffs[key] = modes, h[f"_{store}_lmn"][()]
            profiles = []
            for key in ("_iota", "_pressure"):
                if h[key+"/__class__"][()].decode().split(".")[-1] != "HermiteSplineProfile":
                    raise ValueError("Require saved HermiteSplineProfile")
                x, p = h[key+"/_knots"][()], h[key+"/_params"][()]
                profiles.append(CubicHermiteSpline(x, p[:len(x)], p[len(x):]))
            self.iota, self.pressure = profiles
        if np.any(self.iota(np.linspace(0, 1, 2001)) <= 0):
            raise ValueError("Require positive rotational transform")
        # Exact piecewise-polynomial integral of rho*iota, rather than a
        # tabulated approximate psi label, preserves the saved H5 flux map.
        c = np.zeros((5, len(self.iota.x)-1))
        c[:4] += self.iota.c
        c[1:] += self.iota.c*self.iota.x[:-1]
        self.anti = PPoly(c, self.iota.x).antiderivative()
        self.integral = float(self.anti(1)-self.anti(0))
        curve = self.raw(1., np.linspace(0, TWO_PI, 4097))
        area = abs(trapezoid(curve["R"]*curve["Z_t"], dx=TWO_PI/4096))
        self.a = float(np.sqrt(area/np.pi))
        self.Bn = self.Psi/(np.pi*self.a**2)
        meta = self.path.parent/"metadata.json"
        self.metadata = json.loads(meta.read_text(encoding="utf-8")) if meta.is_file() else {}
        if "equilibrium_sha256" in self.metadata and self.metadata["equilibrium_sha256"] != self.sha:
            raise ValueError("H5 does not match metadata hash")
        if "a_DESC_m" in self.metadata and not np.isclose(self.a, self.metadata["a_DESC_m"], rtol=1e-8):
            raise ValueError("Independent boundary-area normalization disagrees with DESC metadata")

    def u(self, rho):
        return (self.anti(rho)-self.anti(0))/self.integral

    def rho(self, u):
        return float(brentq(lambda r: float(self.u(r))-u, 0., 1., xtol=2e-14))

    def s_derivatives(self, rho):
        s = np.sqrt(self.u(rho))
        up = rho*self.iota(rho)/self.integral
        upp = (self.iota(rho)+rho*self.iota(rho, 1))/self.integral
        return float(s), float(up/(2*s)), float(upp/(2*s)-up**2/(4*s**3))

    def amplitudes(self, rho, key, dr=0):
        modes = self.modes[key]
        ms = np.unique(modes[:, 1])
        values = self.coeffs[key]*np.array([radial(l, m, rho, dr) for l, m, _ in modes])
        return ms, np.array([values[modes[:, 1] == m].sum() for m in ms])

    def raw(self, rho, theta):
        out = {}
        for key in ("R", "Z", "lambda"):
            amp = {dr: self.amplitudes(rho, key, dr) for dr in range(3)}
            for suffix, (dr, dt) in DERIVATIVES.items():
                ms, coeff = amp[dr]
                out[key+suffix] = trig(ms, theta, dt)@coeff
        return out

    def at_pest(self, rho, angle, transform=True):
        ms, coeff = self.amplitudes(rho, "lambda")
        t = invert_angle(lambda t: (trig(ms, t)@coeff, trig(ms, t, 1)@coeff), angle)
        raw = self.raw(rho, t)
        raw["theta"] = t
        return to_pest(raw) if transform else raw
