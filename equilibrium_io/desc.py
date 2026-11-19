"""Axisymmetric DESC 0.17.3 H5 -> COCOS-1 G-EQDSK, without DESC/JAX.

Supports the Fourier-Zernike / HermiteSplineProfile archives used in the
Chance bean-shape study. Unsupported bases/profiles and 3-D inputs fail closed.
The inverse geometry is solved by Newton iteration, not scattered flux fitting.
"""
from pathlib import Path
import warnings

from matplotlib.path import Path as Polygon
import numpy as np
from scipy.constants import mu_0
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
from threadpoolctl import threadpool_limits

from .desc_h5 import H5Equilibrium, radial, trig
from .geqdsk import make_data, fill_exterior, write_geqdsk, padded_bounds


class GeometryInterpolationError(ValueError):
    """The acceleration spline needs more radial knots."""


class SpectralGeometry:
    """Radial cubic interpolation of exact Fourier amplitudes; exact theta modes."""
    def __init__(self, eq, radial_nodes=1025):
        self.eq = eq
        self.splines = {}
        rho = np.linspace(0, 1, radial_nodes)
        self.ms = {}
        for key in ("R", "Z", "lambda"):
            modes = eq.modes[key]
            ms = np.unique(modes[:, 1])
            amplitudes = np.zeros((len(rho), len(ms)))
            lookup = {m: i for i, m in enumerate(ms)}
            for (l, m, _), c in zip(modes, eq.coeffs[key]):
                if c:
                    amplitudes[:, lookup[m]] += c*radial(l, m, rho)
            self.ms[key] = ms
            self.splines[key] = CubicSpline(rho, amplitudes, axis=0)
        self.radial_nodes = radial_nodes
        self.axis = np.array([self.splines[k](0)@trig(self.ms[k], np.array([0.]))[0]
                              for k in ("R", "Z")])
        theta = np.linspace(0, 2*np.pi, 1024, endpoint=False)
        self.boundary = self.exact_tensor(np.array([1.]), theta)
        self.polygon = np.column_stack((self.boundary["R"][0], self.boundary["Z"][0]))
        sample = self.tensor(np.linspace(1/64, 1, 64), theta[::4])
        coords = np.column_stack((sample["R"].ravel(), sample["Z"].ravel()))
        self.seed_rho = np.repeat(np.linspace(1/64, 1, 64), 256)
        self.seed_theta = np.tile(theta[::4], 64)
        self.tree = cKDTree(coords)
        # Independent exact Zernike evaluation checks the acceleration spline.
        error = 0.
        for r in np.linspace(.00731, .99973, 23):
            exact = eq.raw(float(r), theta[::8])
            approx = self.tensor(np.array([r]), theta[::8])
            error = max(error, *(float(np.max(abs(approx[k][0]-exact[k]))) for k in ("R", "Z")))
        self.sampled_interpolation_error_m = error
        if error > 2e-8:
            raise GeometryInterpolationError(f"Radial geometry interpolation error {error:g} m; increase radial_nodes")
        check = self.tensor(np.linspace(.001, 1, 65), theta[::2])
        det = check["R_r"]*check["Z_t"]-check["R_t"]*check["Z_r"]
        if np.any(det >= 0):
            raise ValueError("Require nested right-handed DESC mapping (clockwise theta)")

    def tensor(self, rho, theta):
        out = {}
        for key, spline in self.splines.items():
            a = spline(rho)
            out[key] = a@trig(self.ms[key], theta).T
            out[key+"_r"] = spline(rho, 1)@trig(self.ms[key], theta).T
            out[key+"_t"] = a@trig(self.ms[key], theta, 1).T
        return out

    def paired(self, rho, theta):
        out = {}
        for key in ("R", "Z"):
            spline, ms = self.splines[key], self.ms[key]
            a, basis = spline(rho), trig(ms, theta)
            out[key] = np.sum(a*basis, axis=1)
            out[key+"_r"] = np.sum(spline(rho, 1)*basis, axis=1)
            out[key+"_t"] = np.sum(a*trig(ms, theta, 1), axis=1)
        return out

    def exact_tensor(self, rho, theta):
        """Analytic Zernike values and derivatives for profiles/boundary fields."""
        out = {}
        for key in ("R", "Z", "lambda"):
            ms, modes = self.ms[key], self.eq.modes[key]
            lookup = {m: i for i, m in enumerate(ms)}
            for dr, suffix in ((0, ""), (1, "_r"), (2, "_rr")):
                if key == "lambda" and dr == 2:
                    continue
                a = np.zeros((len(rho), len(ms)))
                for (l, m, _), c in zip(modes, self.eq.coeffs[key]):
                    if c:
                        a[:, lookup[m]] += c*radial(l, m, rho, dr)
                out[key+suffix] = a@trig(ms, theta).T
                if dr < 2:
                    out[key+("_t" if dr == 0 else "_rt")] = a@trig(ms, theta, 1).T
        return out

    def invert(self, points):
        points = np.asarray(points)
        _, index = self.tree.query(points)
        rho, theta = self.seed_rho[index].copy(), self.seed_theta[index].copy()
        maximum_residual = 0.
        for first in range(0, len(points), 4096):
            section = slice(first, first+4096)
            p, r, t = points[section], rho[section], theta[section]
            at_axis = np.linalg.norm(p-self.axis, axis=1) < 1e-11
            active = ~at_axis
            for iteration in range(35):
                if not np.any(active):
                    break
                d = self.paired(r[active], t[active])
                residual = np.column_stack((d["R"], d["Z"]))-p[active]
                error = np.linalg.norm(residual, axis=1)
                ids = np.flatnonzero(active)
                done = error < 2e-11
                active[ids[done]] = False
                if np.all(done):
                    break
                ids = ids[~done]
                det = d["R_r"]*d["Z_t"]-d["R_t"]*d["Z_r"]
                dr = (d["Z_t"]*residual[:, 0]-d["R_t"]*residual[:, 1])/det
                dt = (-d["Z_r"]*residual[:, 0]+d["R_r"]*residual[:, 1])/det
                r[ids] = np.clip(r[ids]-np.clip(dr[~done], -.1, .1), 1e-9, 1.)
                t[ids] = (t[ids]-np.clip(dt[~done], -.5, .5)) % (2*np.pi)
            if np.any(active):
                raise ValueError(f"DESC inverse map did not converge for {active.sum()} interior points")
            r[at_axis] = 0.
            check = self.paired(r, t)
            residual = np.hypot(check["R"]-p[:, 0], check["Z"]-p[:, 1])
            maximum_residual = max(maximum_residual, float(residual.max(initial=0)))
            rho[section], theta[section] = r, t
        return rho, theta, maximum_residual

    def f_samples(self, rho):
        theta = np.linspace(0, 2*np.pi, 512, endpoint=False)
        d = self.exact_tensor(rho, theta)
        det = d["R_r"]*d["Z_t"]-d["R_t"]*d["Z_r"]
        det_r = (d["R_rr"]*d["Z_t"]+d["R_r"]*d["Z_rt"]
                 -d["R_rt"]*d["Z_r"]-d["R_t"]*d["Z_rr"])
        jac = -d["R"]*det
        jac_r = -d["R_r"]*det-d["R"]*det_r
        if np.any(jac <= 0):
            raise ValueError("Invalid DESC Jacobian")
        f = self.eq.Psi*np.asarray(rho)[:, None]/np.pi*d["R"]**2*(1+d["lambda_t"])/jac
        mean = np.sum(f*jac, axis=1)/np.sum(jac, axis=1)
        spread = np.sqrt(np.sum((f-mean[:, None])**2*jac, axis=1)/np.sum(jac, axis=1))
        a = np.sum(d["R"]**2*(1+d["lambda_t"]), axis=1)
        ar = np.sum(2*d["R"]*d["R_r"]*(1+d["lambda_t"])
                    +d["R"]**2*d["lambda_rt"], axis=1)
        denominator, derivative = np.sum(jac, axis=1), np.sum(jac_r, axis=1)
        rho = np.asarray(rho)
        fr = self.eq.Psi/np.pi*((a+rho*ar)/denominator-rho*a*derivative/denominator**2)
        return mean, spread, fr


def desc_to_data(filename, *, nr=257, nz=None, radial_nodes=1025, padding=0.02):
    """Sample an accepted axisymmetric DESC H5 on an R,Z grid.

    F is the flux-surface average of R*Bphi. Its finite angular variation is
    recorded. Pressure and q come directly from the saved profiles. FF' is
    differentiated analytically from the surface average, not inferred from GS.
    """
    if isinstance(nr, bool) or not isinstance(nr, int) or not 16 <= nr <= 999:
        raise ValueError("nr must be an integer in 16..999")
    eq = H5Equilibrium(filename)
    if eq.version != "0.17.3":
        raise ValueError(f"This archive adapter is validated for DESC 0.17.3, received {eq.version}")
    if eq.metadata and not (eq.metadata.get("accepted") is True and eq.metadata.get("quality_passed") is True):
        raise ValueError("DESC metadata does not describe an accepted quality-passed equilibrium")
    with threadpool_limits(limits=1):
        if isinstance(radial_nodes, bool) or not isinstance(radial_nodes, int) or radial_nodes < 65:
            raise ValueError("radial_nodes must be an integer >=65")
        nodes = radial_nodes
        while True:
            try:
                geo = SpectralGeometry(eq, nodes)
                break
            except GeometryInterpolationError:
                if nodes >= 8193:
                    raise
                nodes = min(8193, 2*(nodes-1)+1)
        b = geo.polygon
        lower, upper = padded_bounds(b, padding)
        nz = max(16, int(round(nr*(upper[1]-lower[1])/(upper[0]-lower[0])))) if nz is None else nz
        if isinstance(nz, bool) or not isinstance(nz, int) or not 16 <= nz <= 999:
            raise ValueError("nz must be an integer in 16..999")
        r, z = np.linspace(lower[0], upper[0], nr), np.linspace(lower[1], upper[1], nz)
        rr, zz = np.meshgrid(r, z, indexing="ij")
        points = np.column_stack((rr.ravel(), zz.ravel()))
        inside = Polygon(b).contains_points(points, radius=1e-12)
        rho, theta, inverse_residual = geo.invert(points[inside])
        span = eq.Psi*eq.integral/np.pi
        psi = np.full(len(points), np.nan)
        psi[inside] = eq.u(rho)*span
        d = geo.boundary
        det = d["R_r"][0]*d["Z_t"][0]-d["R_t"][0]*d["Z_r"][0]
        psi_r = eq.Psi/np.pi*float(eq.iota(1.))
        gradient = psi_r*np.column_stack((d["Z_t"][0], -d["R_t"][0]))/det[:, None]
        psi = fill_exterior(points, inside, psi, b, gradient, span)
        u = np.linspace(0., 1., nr)
        profile_rho = np.array([eq.rho(float(v)) for v in u])
        pr = profile_rho.copy(); pr[0] = 1e-3
        f, spread, fr = geo.f_samples(pr)
        tiny_r = np.array([2e-3, 3e-3])
        tiny_f, _, tiny_fr = geo.f_samples(tiny_r)
        # Even-in-rho axis limit, O(rho^6); avoids the 0/0 Jacobian on axis.
        f[0] = 1.5*f[0]-.6*tiny_f[0]+.1*tiny_f[1]
        dp = eq.pressure(profile_rho, 1)
        dpsi = eq.Psi/np.pi*profile_rho*eq.iota(profile_rho)
        ffprime = np.divide(f*fr, dpsi, out=np.zeros_like(f), where=dpsi != 0)
        frr_axis = 1.5*fr[0]/pr[0]-.6*tiny_fr[0]/tiny_r[0]+.1*tiny_fr[1]/tiny_r[1]
        ffprime[0] = f[0]*frr_axis/(eq.Psi/np.pi*float(eq.iota(0.)))
        pprime = np.divide(dp, dpsi, out=np.zeros_like(dp), where=dpsi != 0)
        if abs(dp[0]) > 1e-8*max(1., abs(float(eq.pressure(0.)))):
            raise ValueError("Pressure is not regular on the magnetic axis")
        pprime[0] = float(eq.pressure(0., 2))/(eq.Psi/np.pi*float(eq.iota(0.)))
        jac = -d["R"][0]*det
        factor = eq.Psi/np.pi/jac
        current = float(np.mean(factor*eq.iota(1.)*(d["R_t"][0]**2+d["Z_t"][0]**2))*2*np.pi/mu_0)
        rcentr = float((b[:, 0].min()+b[:, 0].max())/2)
        data = make_data(r, z, psi.reshape(nr, nz), axis=geo.axis, span=span,
            rcentr=rcentr, bcentr=float(f[-1]/rcentr), current=current,
            fpol=f, pres=eq.pressure(profile_rho), ffprim=ffprime, pprime=pprime,
            qpsi=1/eq.iota(profile_rho), boundary=b,
            header="DESC bean COCOS=1 Wb/rad interior-valid")
        metadata = dict(schema=1, source="DESC axisymmetric H5", desc_version=eq.version,
            source_h5_sha256=eq.sha, cocos=1, flux_units="Wb/rad",
            field_convention="BR=+psi_Z/R; BZ=-psi_R/R; Bphi=F/R",
            source_resolution=eq.resolution, grid=[nr, nz],
            grid_padding_fraction=float(padding),
            exterior="normal_linear_continuation_NOT_a_vacuum_solution",
            valid_region="inside supplied LCFS; exterior must not be used for vacuum physics",
            outside_grid_points=int(np.count_nonzero(~inside)),
            radial_geometry_nodes=geo.radial_nodes,
            sampled_radial_interpolation_error_m=geo.sampled_interpolation_error_m,
            inverse_map_max_residual_m=inverse_residual,
            F_surface_average="DESC Jacobian-weighted average of R*Bphi",
            F_surface_relative_std_max=float(np.max(spread/abs(f))),
            Ip_A=current, original_metadata=eq.metadata)
        if "Ip_A" in eq.metadata:
            metadata["Ip_relative_difference_from_DESC"] = current/eq.metadata["Ip_A"]-1
        return data, metadata


def export_desc_gfile(source, destination, *, nr=257, nz=None, radial_nodes=1025,
                      padding=0.02, overwrite=False):
    path = Path(destination)
    if not overwrite and (path.exists() or Path(str(path)+".json").exists()):
        raise FileExistsError(path)
    data, metadata = desc_to_data(source, nr=nr, nz=nz, radial_nodes=radial_nodes, padding=padding)
    warnings.warn("DESC gfile exterior is normal continuation, not a vacuum solution. "
                  "See the .json sidecar.", UserWarning, stacklevel=2)
    return write_geqdsk(data, path, overwrite=overwrite, metadata=metadata)
