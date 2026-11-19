"""Export a solved P3 GS equilibrium as COCOS-1 G-EQDSK."""
from dataclasses import replace
from pathlib import Path
import warnings

import numpy as np

from .ballooning_adapter import build_ballooning_adapter
from .post import (build_post_data, compute_plasma_current, P3FieldEvaluator,
                   evaluate_q_geometry, evaluate_q_from_g)
from equilibrium_io.geqdsk import fill_exterior, write_geqdsk, padded_bounds


def export_gfile(equilibrium, filename, *, nr=257, nz=None, padding=0.02, overwrite=False):
    """Write a gfile plus mandatory provenance sidecar, without changing the solve.

    psi_export = sign(native_span)*(psi_native-native_axis), in Wb/rad.
    F and q retain their values; p', FF' use the positive exported flux span.
    Ip_export = -sign(native_span)*Ip_native, consistent with COCOS 1.
    The exterior is a normal linear continuation, not a physical vacuum field.
    Only converged equilibria are accepted. If q was skipped, it is computed here.
    """
    path = Path(filename)
    if not overwrite and (path.exists() or Path(str(path)+".json").exists()):
        raise FileExistsError(path)
    if not equilibrium.converged:
        raise ValueError("Cannot export an unconverged GS equilibrium")
    if equilibrium.q_result is None:
        geometry = evaluate_q_geometry(equilibrium.mesh, equilibrium.state,
                                       equilibrium.problem.q_options)
        equilibrium = replace(equilibrium, q_result=evaluate_q_from_g(geometry, equilibrium.f_profile))
    data, _, diagnostics = build_ballooning_adapter(equilibrium, nr=nr, nz=nz)
    sign = float(np.sign(equilibrium.state.dpsi))
    # The solver/postprocessor uses piecewise-linear source/F interpolation.
    # The in-memory ballooning adapter's PCHIP is for its own geometry workflow.
    grid = equilibrium.problem.psi_n
    data.fpol = np.interp(data.psin, grid, equilibrium.f_profile.f)
    data.pres = np.interp(data.psin, grid, equilibrium.pressure)
    data.pprime = np.interp(data.psin, grid, equilibrium.pprime)/abs(equilibrium.state.dpsi)
    data.ffprim = np.interp(data.psin, grid, equilibrium.f_profile.ffprime)/abs(equilibrium.state.dpsi)
    native_current = compute_plasma_current(build_post_data(equilibrium))
    data.current = -sign*native_current
    data.header = "P3 GS COCOS=1 Wb/rad interior-valid"
    boundary = np.column_stack((data.rbbbs, data.zbbbs))
    lower, upper = padded_bounds(boundary, padding)
    data.rgrid = np.linspace(lower[0], upper[0], data.nw)
    data.zgrid = np.linspace(lower[1], upper[1], data.nh)
    data.rleft, data.rdim = float(lower[0]), float(upper[0]-lower[0])
    data.zmid, data.zdim = float((lower[1]+upper[1])/2), float(upper[1]-lower[1])
    rr, zz = np.meshgrid(data.rgrid, data.zgrid, indexing="ij")
    points = np.column_stack((rr.ravel(), zz.ravel()))
    evaluator = P3FieldEvaluator(equilibrium.mesh)
    _, gradients, _ = evaluator.evaluate(equilibrium.psi, boundary)
    values = np.empty(len(points))
    inside = np.empty(len(points), dtype=bool)
    for first in range(0, len(points), 32768):
        section = slice(first, first+32768)
        op, valid = evaluator.value_operator(points[section])
        values[section] = sign*(op@equilibrium.psi-equilibrium.state.psi_axis)
        inside[section] = valid
    data.psirz = fill_exterior(points, inside, values, boundary,
                              sign*gradients, data.sibry).reshape(data.nw, data.nh)
    warnings.warn("Gfile exterior uses normal continuation; no vacuum field is supplied. "
                  "See the .json sidecar.", UserWarning, stacklevel=2)
    metadata = dict(schema=1, source="gs_fixed.solve_gs", cocos=1,
        flux_units="Wb/rad", field_convention="BR=+psi_Z/R; BZ=-psi_R/R; Bphi=F/R",
        native_axis=diagnostics.nativeAxis, native_span=diagnostics.nativeSpan,
        native_current_A=float(native_current), native_psi_scale=sign,
        native_poloidal_field_scale=-sign, exported_current_A=float(data.current),
        converged=True, grid=[data.nw, data.nh], grid_padding_fraction=float(padding),
        exterior="normal_linear_continuation_NOT_a_vacuum_solution",
        valid_region="inside supplied LCFS; exterior must not be used for vacuum physics",
        outside_grid_points=int(np.count_nonzero(~inside)),
        axis_flux_offset_Wb_per_rad=float(sign*(equilibrium.q_result.axis.psi
                                               -equilibrium.state.psi_axis)),
        q_edge="extrapolated by existing GS postprocessor")
    return write_geqdsk(data, path, overwrite=overwrite, metadata=metadata)
