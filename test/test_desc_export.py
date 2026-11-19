"""Manufactured Fourier-Zernike geometry: exact flux/profile transformations.

This circular fixture tests field projection and is NOT claimed to solve GS.
"""
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from equilibrium_io.desc import desc_to_data, SpectralGeometry
from equilibrium_io.desc_h5 import H5Equilibrium


def circular_fixture(path):
    rho = np.linspace(0, 1, 17)
    with h5py.File(path, "w") as h:
        for key, value in {"_N": 0, "_L": 2, "_M": 1, "_Psi": 1.2}.items():
            h[key] = value
        h["__version__"] = np.bytes_("0.17.3")
        modes = np.array([[0, 0, 0], [1, 1, 0], [1, -1, 0]])
        for key, values in (("R", [3., .4, 0]), ("Z", [0, 0, -.4]), ("L", [0, 0, 0])):
            h[f"_{key}_basis/__class__"] = np.bytes_("desc.basis.FourierZernikeBasis")
            h[f"_{key}_basis/_modes"] = modes
            h[f"_{key}_lmn"] = values
        for key, values, derivatives in (("iota", np.full_like(rho, .5), np.zeros_like(rho)),
                                         ("pressure", 1000*(1-rho**2), -2000*rho)):
            h[f"_{key}/__class__"] = np.bytes_("desc.profiles.HermiteSplineProfile")
            h[f"_{key}/_knots"] = rho
            h[f"_{key}/_params"] = np.r_[values, derivatives]


class DESCExportTests(unittest.TestCase):
    def test_exact_geometry_profiles_and_flux_factor(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"equilibrium.h5"
            circular_fixture(path)
            d, meta = desc_to_data(path, nr=33, nz=41)
            expected_span = 1.2*.5/(2*np.pi)
            self.assertAlmostEqual(d.sibry, expected_span, places=13)
            np.testing.assert_allclose(d.qpsi, 2., rtol=1e-14)
            np.testing.assert_allclose(d.pres, 1000*(1-d.psin), atol=1e-10)
            np.testing.assert_allclose(d.pprime, -1000/expected_span, rtol=1e-13)
            # Analytic Jacobian-weighted projection of R*Bphi.
            k = 1.2/(np.pi*.4**2)
            expected_f = k*(3+.4**2*d.psin/(2*3))
            np.testing.assert_allclose(d.fpol, expected_f, rtol=2e-12)
            np.testing.assert_allclose(d.ffprim, expected_f*k*.4**2/(6*expected_span), rtol=5e-10)
            rr, zz = np.meshgrid(d.rgrid, d.zgrid, indexing="ij")
            u = ((rr-3)**2+zz**2)/.4**2
            inside = u < .99
            np.testing.assert_allclose(d.psirz[inside], expected_span*u[inside], atol=1e-10)
            self.assertLess(meta["inverse_map_max_residual_m"], 3e-11)
            with h5py.File(path, "a") as h:
                h["_N"][...] = 1
            with self.assertRaises(ValueError):
                desc_to_data(path)


if __name__ == "__main__":
    unittest.main()
