"""File interchange, native-flux signs, exact P3 values and analytic eigenvalues.

Run: python -m unittest test.test_geqdsk
FreeQDSK is an optional independent reader used when installed.
"""
from dataclasses import replace
import gzip
from pathlib import Path
import tempfile
import unittest
import warnings

import numpy as np
from scipy.interpolate import PchipInterpolator

from equilibrium_io.geqdsk import make_data, read_geqdsk, write_geqdsk
from gs_fixed.post import P3FieldEvaluator
from gs_fixed.solver import solve_gs
from test.test_solver import make_problem


class GFileTests(unittest.TestCase):
    def test_asymmetric_roundtrip_and_independent_reader(self):
        r, z = np.linspace(2, 4, 17), np.linspace(-1, 1, 23)
        rr, zz = np.meshgrid(r, z, indexing="ij")
        u = np.linspace(0, 1, len(r))
        data = make_data(r, z, rr**3+2*zz+rr*zz, axis=(3, 0), span=.4,
            rcentr=3, bcentr=-2, current=-1.2e6, fpol=-6+u, pres=1e4*(1-u),
            ffprim=-2-u, pprime=-1e4+u, qpsi=-2-u,
            boundary=np.array([[2.2, -.8], [3.8, -.8], [3.8, .8], [2.2, .8]]))
        with tempfile.TemporaryDirectory() as directory:
            for name in ("gtest", "test.geqdsk.gz"):
                path = Path(directory)/name
                write_geqdsk(data, path)
                result = read_geqdsk(path)
                self.assertEqual(result.psirz.shape, (17, 23))
                for key in ("psirz", "qpsi", "fpol", "pres", "ffprim", "pprime"):
                    np.testing.assert_allclose(getattr(result, key), getattr(data, key), rtol=6e-10, atol=1e-8)
                with self.assertRaises(FileExistsError):
                    write_geqdsk(data, path)
                try:
                    from freeqdsk import geqdsk
                except ImportError:
                    continue
                opener = gzip.open if path.suffix == ".gz" else open
                with opener(path, "rt") as stream:
                    independent = geqdsk.read(stream)
                np.testing.assert_allclose(independent.psi, data.psirz, rtol=6e-10, atol=1e-8)
                np.testing.assert_allclose(independent.ffprime, data.ffprim, rtol=6e-10)
            bad = replace(data, psirz=np.full_like(data.psirz, np.nan))
            with self.assertRaises(ValueError):
                write_geqdsk(bad, Path(directory)/"invalid")
            self.assertFalse((Path(directory)/"invalid").exists())
            with self.assertRaises(ValueError):
                write_geqdsk(replace(data, rdim=99), Path(directory)/"invalid")

    def test_solved_gs_export_and_flux_reversal(self):
        from eq_process import Equilibrium, check_cocos
        problem = make_problem(nr=17, nz=21, contour_nr=96)
        u = problem.psi_n
        eq = solve_gs(problem, -1e6*(1-u**1.5)**2, ffprime=2*(1-u))
        self.assertTrue(eq.converged)
        evaluator = P3FieldEvaluator(eq.mesh)
        with tempfile.TemporaryDirectory() as directory, warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            path = eq.export_gfile(Path(directory)/"gs.geqdsk", nr=65, nz=83)
            d = read_geqdsk(path)
            self.assertEqual(check_cocos(d), 1)
            rr, zz = np.meshgrid(d.rgrid, d.zgrid, indexing="ij")
            values, _, ids = evaluator.evaluate(eq.psi, np.column_stack((rr.ravel(), zz.ravel())))
            inside = ids >= 0
            expected = np.sign(eq.state.dpsi)*(values-eq.state.psi_axis)
            np.testing.assert_allclose(d.psirz.ravel()[inside], expected[inside], atol=6e-10)
            np.testing.assert_allclose(d.pprime, np.interp(d.psin, u, eq.pprime)/abs(eq.state.dpsi), rtol=6e-10, atol=1e-9)
            self.assertTrue(np.all(d.psirz.ravel()[~inside] >= d.sibry-1e-9))
            self.assertTrue(np.all(d.qpsi > 0))
            # The export is usable through the existing gfile-facing API.
            loaded = Equilibrium.from_geqdsk(path)
            surface = loaded.extract_flux_surface(.5)
            self.assertGreater(surface.R.size, 10)
            # Native opposite flux orientation is canonically converted once.
            state = replace(eq.state, psi=-eq.psi, psi_axis=-eq.state.psi_axis,
                            dpsi=-eq.state.dpsi, axis_mode="min", final_rhs=-eq.state.final_rhs)
            q = replace(eq.q_result, axis=replace(eq.q_result.axis, psi=-eq.q_result.axis.psi))
            reversed_eq = replace(eq, state=state, q_result=q)
            rev = read_geqdsk(reversed_eq.export_gfile(Path(directory)/"reverse", nr=65, nz=83))
            np.testing.assert_allclose(rev.psirz, d.psirz, atol=1e-12)
            self.assertAlmostEqual(rev.current/d.current, 1., places=9)
            self.assertEqual(check_cocos(rev), 1)
            with self.assertRaises(ValueError):
                replace(eq, inner_result=replace(eq.inner_result, converged=False)).export_gfile(Path(directory)/"bad")
            self.assertFalse((Path(directory)/"bad").exists())
            noq = replace(eq, q_result=None)
            self.assertTrue(noq.export_gfile(Path(directory)/"computed_q", nr=33).is_file())

    def test_ballooning_eigenvalue_normalization(self):
        from sturm_liouville import Sturm_Liouville
        t = np.linspace(-np.pi/2, np.pi/2, 65)
        result = Sturm_Liouville(t, np.ones_like(t), np.zeros_like(t),
                                np.ones_like(t), numEigenvalues=1)
        self.assertAlmostEqual(float(result.eigval[0]), -1., places=6)


if __name__ == "__main__":
    unittest.main()
