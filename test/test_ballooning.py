"""Analytic finite-difference eigenvalue and coefficient identities."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ballooning.desc import mode, coefficients
from equilibrium_io.desc_h5 import H5Equilibrium
from test.test_desc_export import circular_fixture


class BallooningTests(unittest.TestCase):
    def test_fd_eigenvalue_and_energy(self):
        n = 256
        t = np.linspace(0, np.pi, n+1)
        result, x = mode(np.ones_like(t), np.zeros_like(t), np.ones_like(t), t[1]-t[0])
        exact_discrete = -4*np.sin(np.pi/(2*n))**2/(t[1]-t[0])**2
        self.assertAlmostEqual(result["lambda"], exact_discrete, places=10)
        self.assertLess(abs(result["rayleigh_gap"]), 1e-10)
        self.assertLess(result["matrix_residual_scaled"], 1e-12)
        np.testing.assert_allclose(x[[0, -1]], 0)

    def test_spectral_coefficient_identities(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"equilibrium.h5"
            circular_fixture(path)
            eq = H5Equilibrium(path)
            rho = eq.rho(.5)
            t = np.linspace(-4*np.pi, 4*np.pi, 1025)
            c = coefficients(eq, eq.at_pest(rho, t), rho, t, theta0=.3)
            self.assertTrue(np.all(c["g"] > 0))
            self.assertTrue(np.all(c["f"] > 0))
            self.assertLess(c["metric_identity_error"], 1e-12)
            # This concentric, constant-iota fixture has I=0; its relative
            # identity residual divides roundoff by roundoff. Check I itself.
            self.assertLess(float(np.max(abs(c["I"]))), 1e-12)


if __name__ == "__main__":
    unittest.main()
