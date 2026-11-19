"""Analytic integration, flux conventions, extensions and solved-GS checks.

Run from the project root: python -m test.test_post_global
The manufactured rectangle tests integration, not a fixed-boundary GS PDE.
"""
from dataclasses import replace
import numpy as np
from gs_fixed.fem import MU0, assemble_rhs
from gs_fixed.mesh import P1Mesh, generate_p3_nodes
from gs_fixed.profiles import reconstruct_g_from_ffprime, FFprimeModel, evaluate_ffprime_model
from gs_fixed.solver import (GSEquilibrium, GSState, PicardResult, prepare_gs_problem,
                             solve_gs, PicardOptions, QOptions, PQOptions, ParallelOptions)
from gs_fixed.profiles import QTarget
from gs_fixed.post import (build_post_data, compute_global_parameters, evaluate_magnetic_field,
                           compute_plasma_current, compute_beta_pol, compute_li3, compute_beta_n)


def manufactured_equilibrium():
    mesh = generate_p3_nodes(P1Mesh(np.array([[2.,-.5], [4.,-.5], [4.,.5], [2.,.5]]),
                                     np.array([[0,1,2], [0,2,3]])))
    problem = prepare_gs_problem(mesh=mesh, r0=3., b0=3., a_minor=1., n_profile=101)
    r, z = mesh.points.T
    psi = 2. - (r-3.)**2 - z**2
    pp = np.full_like(problem.psi_n, -2.e4)
    ff = np.full_like(problem.psi_n, .05)
    rhs = assemble_rhs(mesh, psi, problem.quad, pp, ff, -2.)
    state = GSState(psi, (psi-2.)/-2., 2., -2., "max", int(np.argmax(psi)),
                    np.array([3.,0.]), rhs, np.zeros(len(problem.bc.free_nodes)), True, 0)
    inner = PicardResult(True, 0, *[np.empty(0) for _ in range(6)], 0., 0.)
    profile = reconstruct_g_from_ffprime(problem.psi_n, ff, 9.)
    return GSEquilibrium(1, problem, state, pp, profile,
                         3.e4 - 2.e4*problem.psi_n, inner)


def expect_value_error(function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except ValueError:
        return
    raise AssertionError("Expected ValueError")


def test_analytic_integrals():
    eq = manufactured_equilibrium()
    default = build_post_data(eq)
    np.testing.assert_allclose(compute_plasma_current(default), default.ip_from_rhs, rtol=2e-14)
    # Closed-form moments on R in [2,4], Z in [-.5,.5].
    expected_volume = 12*np.pi
    expected_pressure = 2*np.pi*(6*3.e4 + 1.25*(-2.e4))
    expected_bp2 = 8*np.pi*(-6 + (109/12)*np.log(2))
    expected_ip = 6*(-2.e4)/(-2) + .05/(-2*MU0)*np.log(2)
    errors = []
    for degree in (6,10,20):
        data = build_post_data(eq, quadrature_degree=degree, element_chunk_size=1)
        errors.append(abs(data.integrate_volume(data.bp_squared)/expected_bp2 - 1))
    assert errors[-1] < errors[0] and errors[-1] < 5e-13
    data = build_post_data(eq, quadrature_degree=20)
    result = compute_global_parameters(data)
    expected_beta_tor = 2*MU0*(expected_pressure/expected_volume)/9
    expected = dict(volume=expected_volume, pressure_integral=expected_pressure,
                    pressure_volume_average=expected_pressure/expected_volume,
                    bp2_volume_integral=expected_bp2, ip=expected_ip,
                    poloidal_magnetic_energy=expected_bp2/(2*MU0),
                    beta_tor=expected_beta_tor,
                    beta_pol=4*expected_pressure/(3*MU0*expected_ip**2),
                    li3=2*expected_bp2/(3*MU0**2*expected_ip**2),
                    beta_n=100*expected_beta_tor*3/(abs(expected_ip)/1.e6),
                    b_axis=eq.f_profile.f[0]/3)
    for name, value in expected.items():
        np.testing.assert_allclose(result[name], value, rtol=6e-13, atol=1e-12, err_msg=name)
    np.testing.assert_allclose(data.integrate_area(), 2., rtol=1e-14)
    np.testing.assert_allclose(data.integrate_volume(2.), 2*expected_volume, rtol=1e-14)
    np.testing.assert_allclose(data.volume_average(1.), 1., rtol=1e-14)
    assert not data.pressure.flags.writeable
    # A new physical definition needs only a function and an optional dict key.
    def compute_mean_b_squared(samples):
        return samples.volume_average(samples.b_squared)
    result["mean_b_squared"] = compute_mean_b_squared(data)
    assert result["mean_b_squared"] > data.volume_average(data.bp_squared)
    expect_value_error(data.integrate_volume, np.ones((2,1)))
    expect_value_error(data.integrate_volume, np.nan)
    expect_value_error(data.integrate_volume, 1j)
    print("Analytic integrals passed; Bp^2 errors at degrees 6/10/20:", errors)


def test_fields_and_conventions():
    eq = manufactured_equilibrium()
    data = build_post_data(eq)
    reversed_state = replace(eq.state, psi=-eq.psi, psi_axis=-2., dpsi=2.,
                             axis_mode="min", final_rhs=-eq.state.final_rhs)
    reverse = replace(eq, state=reversed_state)
    reverse_data = build_post_data(reverse)
    np.testing.assert_allclose(reverse_data.psi_n, data.psi_n, rtol=0, atol=0)
    np.testing.assert_allclose(reverse_data.b_r, -data.b_r, rtol=0, atol=0)
    np.testing.assert_allclose(reverse_data.b_z, -data.b_z, rtol=0, atol=0)
    np.testing.assert_allclose(reverse_data.b_phi, data.b_phi, rtol=0, atol=0)
    np.testing.assert_allclose(compute_plasma_current(reverse_data), -compute_plasma_current(data))
    a, b = compute_global_parameters(data), compute_global_parameters(reverse_data)
    for name in ("beta_tor", "beta_pol", "beta_n", "li3", "volume", "poloidal_magnetic_energy"):
        np.testing.assert_allclose(a[name], b[name], rtol=0, atol=0)
    shifted_state = replace(eq.state, psi=eq.psi+100., psi_axis=102.)
    shifted = build_post_data(replace(eq, state=shifted_state))
    np.testing.assert_allclose(shifted.bp_squared, data.bp_squared, rtol=3e-12, atol=1e-13)
    np.testing.assert_allclose(shifted.pressure, data.pressure, rtol=1e-13)
    r, z = np.broadcast_arrays(np.array([[2.3],[3.2]]), np.array([-.2,.3]))
    field = evaluate_magnetic_field(eq, r, z)
    assert field["psi"].shape == (2,2) and np.all(field["inside"])
    np.testing.assert_allclose(field["psi"], 2-(r-3)**2-z**2, atol=2e-14)
    np.testing.assert_allclose(field["b_r"], 2*z/r, atol=2e-14)
    np.testing.assert_allclose(field["b_z"], -2*(r-3)/r, atol=2e-14)
    np.testing.assert_allclose(field["b_phi"], field["f"]/r, rtol=1e-15)
    outside = evaluate_magnetic_field(eq, [2.5,10.], 0.)
    np.testing.assert_array_equal(outside["inside"], [True,False])
    assert all(np.isnan(value[1]) for name,value in outside.items() if name != "inside")
    scalar = evaluate_magnetic_field(eq, 3., 0.)
    assert scalar["b"].shape == ()
    np.testing.assert_allclose(scalar["bp"], 0., atol=2e-14)
    failures = [dict(element_chunk_size=0), dict(flux_range_tolerance=-1)]
    for options in failures:
        expect_value_error(build_post_data, eq, **options)
    # The existing quadrature API raises TypeError for a noninteger degree.
    try:
        build_post_data(eq, quadrature_degree=True)
    except TypeError:
        pass
    else:
        raise AssertionError("Boolean quadrature degree was accepted")
    bad = replace(eq, pressure=-np.ones_like(eq.pressure))
    expect_value_error(build_post_data, bad)
    failed = replace(eq, inner_result=replace(eq.inner_result, converged=False))
    expect_value_error(build_post_data, failed)
    assert not build_post_data(failed, allow_unconverged=True).converged
    out_of_range = replace(eq, state=replace(eq.state, psi_axis=10., dpsi=-2.))
    expect_value_error(build_post_data, out_of_range)
    # Zero current leaves volume/energy samples usable; current-normalized
    # quantities explicitly reject their undefined normalization.
    zero = replace(data, j_phi=np.zeros_like(data.j_phi), ip_from_rhs=0.)
    assert compute_plasma_current(zero) == 0
    for function in (compute_beta_pol, compute_li3, compute_beta_n):
        expect_value_error(function, zero)
    print("Native flux sign/offset, spatial fields, extensions and invalid inputs passed.")


def test_solved_equilibrium():
    theta = np.arange(80)*2*np.pi/80
    boundary = np.column_stack((3+.8*np.cos(theta+np.arcsin(.25)*np.sin(theta)),
                                1.2*np.sin(theta)))
    problem = prepare_gs_problem(boundary, r0=3., b0=3., a_minor=.8,
                                 nr=17, nz=21, n_profile=101,
                                 q_options=QOptions(contour_nr=80),
                                 picard=PicardOptions(max_iterations=150,
                                     update_tolerance=1e-9, residual_tolerance=1e-9))
    pp = -1.e6*(1-problem.psi_n**1.5)**2
    ff = evaluate_ffprime_model(FFprimeModel(problem.psi_n, np.array([2.,.2,-.1]), 9.)).ffprime
    eq = solve_gs(problem, pp, ffprime=ff)
    result = compute_global_parameters(build_post_data(eq))
    assert result["converged"] and result["ip"] > 0
    assert result["volume"] > 0 and result["pressure_volume_average"] > 0
    assert result["beta_tor"] > 0 and result["beta_pol"] > 0 and result["li3"] > 0
    assert result["ip_rhs_relative_difference"] < 3e-14
    print("Solved GS postprocessing passed:", result)
    qeq = solve_gs(problem, pp, q_target=QTarget(eq.q_result.psi_n, eq.q_result.q),
                    ffprime_model=FFprimeModel(problem.psi_n, np.array([1.,0.,0.]), 9.),
                    outer_options=PQOptions(max_iterations=15, q_relative_tolerance=1e-4,
                        q_rms_log_tolerance=1e-4, regularization=0., verbose=False),
                    parallel=ParallelOptions(max_workers=1))
    qresult = compute_global_parameters(build_post_data(qeq))
    assert qeq.mode == 3 and qresult["converged"]
    assert qresult["ip_rhs_relative_difference"] < 3e-14
    assert np.all(np.isfinite(list(qresult.values())))
    print("Mode-3 equilibrium uses the same postprocessing interface and passes RHS consistency.")


if __name__ == "__main__":
    test_analytic_integrals()
    test_fields_and_conventions()
    test_solved_equilibrium()
