"""Native spectral DESC ballooning coefficients and largest-eigenvalue solver.

Preserves the validated Chance-study arithmetic; no coarse remapping is applied.
"""
import numpy as np
from scipy.constants import mu_0
from scipy.linalg import eigh_tridiagonal
from threadpoolctl import threadpool_limits
from equilibrium_io.desc_h5 import H5Equilibrium

def coefficients(eq, data, rho, theta, theta0=0.):
    """Chance Eq.(3), zeta=q*theta_PEST, native DESC lambda normalization.

    data is expressed in PEST coordinates: lambda=0, derivatives at fixed rho.
    Curvature is evaluated both geometrically and through force balance. A
    coarsened geometry need not satisfy force balance: these routes must NOT
    be silently identified. All vector components use orthonormal (R,phi,Z).
    """
    r = data["R"]
    rr, rt, rrr, rrt, rtt = (data["R"+s] for s in ("_r", "_t", "_rr", "_rt", "_tt"))
    zr, zt, zrr, zrt, ztt = (data["Z"+s] for s in ("_r", "_t", "_rr", "_rt", "_tt"))
    det = rr*zt-rt*zr
    jac = -r*det
    if not np.all(np.isfinite(jac) & (jac > 0)):
        raise ValueError("Nonpositive reconstructed Jacobian")
    jr = -rr*det-r*(rrr*zt+rr*zrt-rrt*zr-rt*zrr)
    jt = -rt*det-r*(rrt*zt+rr*ztt-rtt*zr-rt*zrt)
    zero = np.zeros_like(r)
    grad_rho = np.column_stack((zt, zero, -rt))/det[:, None]
    grad_theta = np.column_stack((-zr, zero, rr))/det[:, None]
    grad_phi = np.column_stack((zero, 1/r, zero))
    io, ir, pr = (float(eq.iota(rho)), float(eq.iota(rho, 1)), float(eq.pressure(rho, 1)))
    tor_r = eq.Psi*rho/np.pi
    psi_r = io*tor_r
    factor = tor_r/jac
    fr, ft = factor*(1/rho-jr/jac), -factor*jt/jac
    bz, bt = factor, factor*io
    B = np.column_stack((bt*rt, bz*r, bt*zt))
    Br = np.column_stack(((fr*io+factor*ir)*rt+bt*rrt, fr*r+factor*rr,
                         (fr*io+factor*ir)*zt+bt*zrt))
    Bt = np.column_stack((io*(ft*rt+factor*rtt), ft*r+factor*rt,
                         io*(ft*zt+factor*ztt)))
    dot = lambda a, b: np.sum(a*b, axis=1)
    b2 = dot(B, B)
    b = np.sqrt(b2)
    unit = B/b[:, None]
    unit_t = Bt/b[:, None]-B*(dot(B, Bt)/b**3)[:, None]
    kappa = (bt/b)[:, None]*unit_t+(bz/b)[:, None]*np.column_stack((-unit[:, 1], unit[:, 0], zero))
    gradb2 = (2*dot(B, Br))[:, None]*grad_rho+(2*dot(B, Bt))[:, None]*grad_theta
    gradp = pr*grad_rho
    total = mu_0*gradp+.5*gradb2
    kappa_f = (total-unit*dot(unit, total)[:, None])/b2[:, None]
    zeta = theta/io
    alpha_d = (-ir*zeta+ir/io*theta0)[:, None]*grad_rho+grad_theta-io*grad_phi
    grad_alpha = -alpha_d/io
    grad_psi = psi_r*grad_rho
    gp2 = dot(grad_psi, grad_psi)
    I = -dot(grad_alpha, grad_psi)/gp2
    metric = 1/gp2+gp2*I**2/b2
    scale = eq.a**3*eq.Bn*rho**2*io**2
    out = {"g": scale*metric*bz, "f": scale*eq.Bn**2/eq.a**2*metric/bz,
           "I": I, "B": B, "Bzeta": bz, "jac": jac, "grad_psi2": gp2,
           "R": data["R"], "Z": data["Z"]}
    for key, curvature in (("geometry", kappa), ("force", kappa_f)):
        out[key+"_normal"] = scale*2*mu_0*pr/psi_r*dot(curvature, grad_psi)/(bz*gp2)
        out[key+"_shear"] = -scale*2*mu_0*pr/psi_r*I*dot(curvature, np.cross(B, grad_psi))/(bz*b2)
        out[key+"_c"] = out[key+"_normal"]+out[key+"_shear"]
        out[key+"_kappa_n"] = dot(curvature, grad_psi)/np.sqrt(gp2)
    # In PEST coordinates I = q'_psi*(theta-theta0)
    #   - q*(R_t R_r+Z_t Z_r)/(psi_r*(R_t^2+Z_t^2)).
    # Differentiate this expression analytically. No FFT truncation or numerical
    # difference of large secular shear terms is needed to obtain local S.
    qp = -ir/io**2/psi_r
    numerator, denominator = rt*rr+zt*zr, rt**2+zt**2
    nt = rtt*rr+rt*rrt+ztt*zr+zt*zrt
    dt = 2*(rt*rtt+zt*ztt)
    iper = -numerator/(io*psi_r*denominator)
    out["S"] = io*bz*(qp-(nt*denominator-numerator*dt)/(io*psi_r*denominator**2))
    out["shear_identity_error"] = float(np.linalg.norm(I-qp*(theta-theta0)-iper)/max(np.linalg.norm(I), 1e-300))
    out["F_Tm"] = B[:, 1]*r
    metric_direct = dot(grad_alpha, grad_alpha)/b2
    out["metric_identity_error"] = float(np.linalg.norm(metric-metric_direct)/np.linalg.norm(metric_direct))
    out["force_balance_residual_relative"] = float(np.linalg.norm((kappa-kappa_f)*b2[:, None])/max(np.linalg.norm(mu_0*gradp), 1e-300))
    if not np.all(np.isfinite(out["g"]) & (out["g"] > 0) & np.isfinite(out["f"]) & (out["f"] > 0)):
        raise ValueError("Invalid ballooning tension/inertia")
    return out


def mode(g, c, f, step, normal=None, shear=None):
    """Same arithmetic face-g FD stencil as DESC, zero Dirichlet endpoints."""
    face = .5*(g[1:]+g[:-1])
    diagonal = (c[1:-1]-(face[1:]+face[:-1])/step**2)/f[1:-1]
    off = face[1:-1]/step**2/np.sqrt(f[1:-2]*f[2:-1])
    n = len(diagonal)
    val, vec = eigh_tridiagonal(diagonal, off, select="i", select_range=(n-1, n-1))
    x = np.r_[0., vec[:, 0]/np.sqrt(f[1:-1]*step), 0.]
    if x[np.argmax(abs(x))] < 0:
        x = -x
    result = energy(g, c, f, step, x, normal, shear)
    result["lambda"] = float(val[0])
    result["rayleigh_gap"] = result["rayleigh"]-float(val[0])
    residual = diagonal*vec[:, 0]-val[0]*vec[:, 0]
    residual[:-1] += off*vec[1:, 0]
    residual[1:] += off*vec[:-1, 0]
    result["matrix_residual_scaled"] = float(np.linalg.norm(residual)/(np.linalg.norm(diagonal*vec[:, 0])+2*np.linalg.norm(off)+abs(val[0])))
    return result, x


def energy(g, c, f, step, x, normal=None, shear=None):
    mass = float(np.sum(f[1:-1]*x[1:-1]**2)*step)
    tension = -float(np.sum(.5*(g[1:]+g[:-1])*np.diff(x)**2)/step)/mass
    drive = float(np.sum(c[1:-1]*x[1:-1]**2)*step)/mass
    out = {"rayleigh": tension+drive, "tension": tension, "drive": drive}
    if normal is not None:
        out.update(normal=float(np.sum(normal[1:-1]*x[1:-1]**2)*step)/mass,
                   shear=float(np.sum(shear[1:-1]*x[1:-1]**2)*step)/mass)
    return out


def solve_desc_ballooning(equilibrium, psi_n, *, theta0=0., theta_max_pi=16,
                          cells_per_turn=512):
    """One surface/phase, geometric curvature, positive lambda = instability.

    This is one sampled eigenvalue; a nonpositive value does not certify global
    stability. Increase domain/resolution and scan surfaces/phases as required.
    theta0 is in PEST radians; the eigenproblem coordinate is zeta=q*theta_PEST.
    """
    if not np.isfinite(psi_n) or not 0 < psi_n < 1:
        raise ValueError("psi_n must be in (0,1)")
    if not np.isfinite(theta0):
        raise ValueError("theta0 must be finite")
    for name, value in (("theta_max_pi", theta_max_pi), ("cells_per_turn", cells_per_turn)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 2:
            raise ValueError(f"{name} must be an integer >=2")
    eq = equilibrium if isinstance(equilibrium, H5Equilibrium) else H5Equilibrium(equilibrium)
    rho = eq.rho(float(psi_n))
    theta = np.linspace(-theta_max_pi*np.pi, theta_max_pi*np.pi,
                        theta_max_pi*cells_per_turn+1)
    with threadpool_limits(limits=1):
        data = eq.at_pest(rho, theta)
        coeff = coefficients(eq, data, rho, theta, theta0)
        step = float((theta[1]-theta[0])/eq.iota(rho))
        result, eigenfunction = mode(coeff["g"], coeff["geometry_c"], coeff["f"], step,
                                     coeff["geometry_normal"], coeff["geometry_shear"])
    return dict(result, psi_n=float(psi_n), theta0=float(theta0),
                theta=theta, eigenfunction=eigenfunction, coefficients=coeff)
