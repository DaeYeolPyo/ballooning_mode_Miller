"""Analytic P3 field/axis/q and contour-resolution verification."""
import numpy as np
from dataclasses import replace
from gs_fixed.solver import prepare_gs_problem,GSState,PicardOptions
from gs_fixed.post import (
    QOptions,P3FieldEvaluator,build_q_contour_cache,evaluate_q_geometry,evaluate_q_from_g,
)
from gs_fixed.profiles import reconstruct_g_from_ffprime


def test_post():
    radius=.8; major=3.
    theta=np.arange(100)*2*np.pi/100
    boundary=np.column_stack((major+radius*np.cos(theta),radius*np.sin(theta)))
    problem=prepare_gs_problem(boundary,r0=major,b0=3.,a_minor=radius,nr=19,nz=19,
                                n_profile=101,q_options=QOptions(contour_nr=80))
    mesh=problem.mesh; R,Z=mesh.points.T
    psi=.5*(1-((R-major)**2+Z**2)/radius**2)
    node=int(np.argmax(psi))
    state=GSState(psi,(psi-.5)/(-.5),.5,-.5,"max",node,mesh.points[node],
                   np.zeros_like(psi),np.zeros(len(problem.bc.free_nodes)),True,0)
    evaluator=P3FieldEvaluator(mesh)
    from gs_fixed.mesh import P1Mesh,generate_p3_nodes
    patch_mesh=generate_p3_nodes(P1Mesh(np.array([[2.,0.],[3.,0.],[3.,1.],[2.,1.]]),
                                     np.array([[0,1,2],[0,2,3]])))
    locator=P3FieldEvaluator(patch_mesh)
    ids,_,_=locator.locate(np.array([[2.25,.25],[2.75,.75],[2.,0.],[3.,1.]]))
    np.testing.assert_array_equal(ids,[1,1,1,1])
    errors=[]
    for nr in (80,160):
        options=QOptions(contour_nr=nr,interior_psi_n=(.1,.3,.5,.7,.9))
        cache=build_q_contour_cache(mesh,options,evaluator)
        geometry=evaluate_q_geometry(mesh,state,options,evaluator=evaluator,cache=cache)
        np.testing.assert_allclose(geometry.axis.point,[major,0.],atol=1e-13)
        np.testing.assert_allclose(geometry.axis.hessian,-np.eye(2)/radius**2,atol=2e-13)
        np.testing.assert_allclose(geometry.q_per_f[0],radius**2/major,atol=1e-14)
        expected=radius**2/np.sqrt(major**2-radius**2*np.array(options.interior_psi_n))
        error=np.max(abs(geometry.q_per_f[1:-1]/expected-1))
        errors.append(error)
        negative=replace(state,psi=-psi,psi_axis=-.5,dpsi=.5,axis_mode="min")
        reversed_geometry=evaluate_q_geometry(mesh,negative,options,evaluator=evaluator,cache=cache)
        np.testing.assert_allclose(reversed_geometry.q_per_f[:-1],geometry.q_per_f[:-1],rtol=1e-12)
        profile=reconstruct_g_from_ffprime(problem.psi_n,np.zeros_like(problem.psi_n),-9.)
        result=evaluate_q_from_g(geometry,profile)
        assert np.all(result.q<0)
        np.testing.assert_allclose(result.q[:-1],-9*geometry.q_per_f[:-1],rtol=1e-14)
        try:
            evaluate_q_geometry(mesh,state,options,evaluator=evaluator,
                                cache=replace(cache,mesh_points=cache.mesh_points+.1))
        except ValueError: pass
        else: raise AssertionError("Wrong-mesh contour cache accepted.")
    assert errors[1]<errors[0] and errors[1]<3e-4
    print("Analytic axis and q-contour checks passed; relative q/F errors:",errors)


if __name__=="__main__":
    test_post()
