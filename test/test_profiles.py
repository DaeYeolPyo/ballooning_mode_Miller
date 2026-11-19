"""Run from the Python root: python -m test.test_profiles."""
import numpy as np
from numpy.polynomial.legendre import legvander
from gs_fixed.profiles import (
    FFprimeModel, FFprimeOptions, FFprimeInitial, QTarget,
    validate_profile_grid, construct_ffprime_basis,
    reconstruct_g_from_ffprime, reconstruct_pressure,
    evaluate_ffprime_model, initialize_ffprime_model,
    evaluate_q_target, normalize_q_target,
)


def test_profiles():
    s = np.linspace(0.0, 1.0, 101)
    dims = dict(r0=3.0, b0=3.0, a_minor=0.8)
    pprime = -1.0e6 * (1 - s**1.5)**2
    q = 1 + 2*s
    for n in range(1, 13):
        expected = legvander(2*s-1, n-1)
        np.testing.assert_allclose(construct_ffprime_basis(s,n,False),expected,atol=3e-15,rtol=2e-14)
        edge = construct_ffprime_basis(s,n)
        np.testing.assert_allclose(edge,(1-s)[:,None]*expected,atol=3e-15,rtol=2e-14)
        np.testing.assert_array_equal(edge[-1],np.zeros(n))

    # Integral of the linear normalized derivative is analytic.
    original_grid = s.copy()
    for fb in (9.0,-9.0):
        profile = reconstruct_g_from_ffprime(s,2*(1-s),fb)
        exact = 81 - 2*(1-s)**2
        np.testing.assert_allclose(profile.g,exact,rtol=3e-15,atol=3e-14)
        np.testing.assert_allclose(profile.f**2,profile.g,rtol=3e-15)
        assert profile.f[-1] == fb and profile.g[-1] == fb**2
        assert profile.sign_f == np.sign(fb) and profile.boundary_error == 0
        assert not profile.ffprime.flags.writeable
        np.testing.assert_array_equal(profile.dg_dpsi_n,4*(1-s))
    np.testing.assert_array_equal(s,original_grid)

    # Integration remains second-order for nonlinear samples, as in MATLAB.
    errors=[]
    for n in (51,101,201):
        x=np.linspace(0,1,n)
        result=reconstruct_g_from_ffprime(x,(1-x)**2,9)
        analytic=81-(2/3)*(1-x)**3
        errors.append(np.max(abs(result.g-analytic)))
    assert errors[0]/errors[1] > 3.9 and errors[1]/errors[2] > 3.9
    nonuniform=np.array([0.,.03,.17,.4,.7,1.])
    pressure=reconstruct_pressure(nonuniform,-2e4*(1-nonuniform),p_edge=5.)
    np.testing.assert_allclose(pressure,5+1e4*(1-nonuniform)**2,rtol=3e-15,atol=2e-11)
    assert pressure[-1] == 5.

    # Independent trial models preserve the basis dimension and original data.
    coefficients=np.array([2.,0.,0.])
    model=FFprimeModel(s,coefficients,9.)
    coefficients[0]=99
    trial=model.with_coefficients([3.,.1,-.1])
    assert model.coefficients[0] == 2 and trial.n_basis == model.n_basis == 3
    np.testing.assert_allclose(evaluate_ffprime_model(model).ffprime,2*(1-s),atol=2e-15)

    # Exact in-basis projection; a separate coarse sample table uses PCHIP.
    supplied=initialize_ffprime_model(
        s,pprime,q,**dims,
        initial=FFprimeInitial(kind="ffprime",psi_n=s,ffprime=2*(1-s)),
    )
    np.testing.assert_allclose(evaluate_ffprime_model(supplied).ffprime,2*(1-s),rtol=0,atol=3e-12)
    short=np.array([.1,.25,.6,.85])
    extrapolated=initialize_ffprime_model(
        s,pprime,q,**dims,
        initial=FFprimeInitial(kind="ffprime",psi_n=short,ffprime=2*(1-short)),
    )
    np.testing.assert_allclose(extrapolated.initialization.source_ffprime,2*(1-s),atol=2e-15)

    for name in ("constant","zero","constant_g"):
        constant=initialize_ffprime_model(s,pprime,**dims,initial=name)
        result=evaluate_ffprime_model(constant)
        np.testing.assert_array_equal(result.ffprime,np.zeros_like(s))
        np.testing.assert_array_equal(result.g,np.full_like(s,81.))
    negative=initialize_ffprime_model(s,pprime,**dims,initial="constant",f_boundary=-9,f_sign=-1)
    assert evaluate_ffprime_model(negative).f[-1] == -9

    # All accepted q representations agree for a linear target.
    targets=(QTarget(short,1+2*short),{"psi_n":short,"q":1+2*short},
             lambda x:1+2*x, q)
    for target in targets:
        np.testing.assert_allclose(evaluate_q_target(target,s),q,atol=3e-15)
    np.testing.assert_array_equal(evaluate_q_target(2.,s),np.full_like(s,2.))
    canonical=normalize_q_target(q,s)
    query=np.array([0.,.123,.77,1.])
    np.testing.assert_allclose(evaluate_q_target(canonical,query),1+2*query,atol=2e-15)

    # Cylindrical initializer: signs, clipping, nonuniform grids and options.
    cylindrical=initialize_ffprime_model(s,pprime,q,**dims)
    negative_cylinder=initialize_ffprime_model(s,pprime,-q,**dims,f_boundary=-9,f_sign=-1)
    np.testing.assert_allclose(cylindrical.coefficients,negative_cylinder.coefficients,rtol=0,atol=0)
    history=cylindrical.initialization.cylindrical
    assert history.final_cylinder_ffprime.shape == s.shape
    assert cylindrical.initialization.source_ffprime.shape == s.shape
    assert len(history.coefficient_history)==6
    assert np.all(history.flux_span_magnitude_history>0)
    np.testing.assert_array_equal(cylindrical.coefficients[1:],np.zeros(2))
    assert evaluate_ffprime_model(cylindrical).minimum_g>0
    custom=initialize_ffprime_model(
        nonuniform,-1e6*(1-nonuniform**1.5)**2,QTarget(short,1+2*short),
        **dims,options=FFprimeOptions(n_basis=4,cylinder_iterations=3,
            cylinder_relaxation=.7,cylinder_maximum_g_drop=.4),
    )
    assert custom.n_basis==4 and len(custom.initialization.cylindrical.coefficient_history)==3
    for amplitude,expected_sign in ((-1e12,1),(1e12,-1)):
        clipped=initialize_ffprime_model(s,amplitude*np.ones_like(s),q,**dims)
        c=clipped.coefficients[0]
        assert np.sign(c)==expected_sign and -81<=c<=.8*81
        assert evaluate_ffprime_model(clipped).minimum_g>0
    minimum_grid=initialize_ffprime_model(
        [0.,1.],[-1e4,0.],1.5,**dims,
    )
    assert evaluate_ffprime_model(minimum_grid).minimum_g>0

    # Invalid grids, inputs, physical branches and inconsistent signs.
    def fails(function,*args,**kwargs):
        try: function(*args,**kwargs)
        except (ValueError,TypeError): return
        raise AssertionError("Invalid input was accepted.")
    for invalid in ([.1,.5,1.],[0.,.5,.5,1.],[0.,1.,np.nan],[0.]):
        fails(validate_profile_grid,invalid)
    fails(validate_profile_grid,nonuniform,require_uniform=True)
    fails(construct_ffprime_basis,s,13)
    fails(construct_ffprime_basis,[-.1,.5],3)
    fails(construct_ffprime_basis,s,True)
    fails(construct_ffprime_basis,s,3,enforce_zero_edge="false")
    fails(reconstruct_g_from_ffprime,s,np.full_like(s,50.),9.)
    fails(reconstruct_g_from_ffprime,s,2*(1-s),0.)
    fails(reconstruct_g_from_ffprime,s,1j*np.ones_like(s),9.)
    fails(reconstruct_g_from_ffprime,s,np.ones(5),9.)
    fails(model.with_coefficients,[1.,2.])
    fails(FFprimeOptions,cylinder_iterations=0)
    fails(FFprimeOptions,cylinder_relaxation=1.1)
    fails(FFprimeOptions,cylinder_maximum_g_drop=1.)
    fails(FFprimeInitial,kind="ffprime")
    fails(QTarget,short,np.zeros_like(short))
    fails(initialize_ffprime_model,s,pprime,np.zeros_like(s),**dims)
    fails(initialize_ffprime_model,s,pprime,q,**dims,f_boundary=-9,f_sign=1)
    fails(initialize_ffprime_model,s,pprime,q,**dims,normalization="free_fb")
    fails(initialize_ffprime_model,s,pprime,q,**dims,
          initial=FFprimeInitial(kind="ffprime",ffprime=[1.,0.]))

    # Existing FEM consumes normalized derivatives without another conversion.
    from gs_fixed.mesh import P1Mesh,generate_p3_nodes
    from gs_fixed.fem import precompute_quadrature,assemble_rhs
    mesh=generate_p3_nodes(P1Mesh(
        np.array([[2.,0.],[3.,0.],[3.,1.],[2.,1.]]),
        np.array([[0,1,2],[0,2,3]]),
    ))
    psi=np.zeros(mesh.n_dof)
    psi[mesh.free_nodes]=.5
    quadrature=precompute_quadrature(10)
    profile=evaluate_ffprime_model(model)
    actual=assemble_rhs(mesh,psi,quadrature,pprime,profile.ffprime,-.5)
    expected=assemble_rhs(mesh,psi,quadrature,pprime,2*(1-s),-.5)
    np.testing.assert_allclose(actual,expected,rtol=3e-15,atol=2e-16)
    print("Profile reconstruction, basis, initialization, q target and FEM integration checks passed.")
    print("Nonlinear-profile trapezoidal refinement errors:",errors)
    print("Cylindrical leading coefficient:",cylindrical.coefficients[0])


if __name__ == "__main__":
    test_profiles()
