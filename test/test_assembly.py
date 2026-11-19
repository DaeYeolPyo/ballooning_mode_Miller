import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, factorized
from unittest.mock import patch
from dataclasses import replace
from gs_fixed.mesh import P1Mesh, generate_p3_nodes, build_p1_mesh
from gs_fixed.fem import (
    MU0, P3_REFERENCE_NODES, precompute_quadrature,
    construct_k_matrix, assemble_k_matrix,
    construct_rhs, assemble_rhs, apply_bc,
)
import gs_fixed.fem as fem



def test_assembly():
    p=np.array([[2.,0.],[3.,0.],[3.,1.],[2.,1.]])
    mesh=generate_p3_nodes(P1Mesh(p,np.array([[0,1,2],[0,2,3]])))
    quad=precompute_quadrature(10)
    K=assemble_k_matrix(mesh,quad)
    assert sparse.isspmatrix_csr(K)
    np.testing.assert_allclose(K.toarray(),K.toarray().T,atol=2e-15)
    np.testing.assert_allclose(K@np.ones(mesh.n_dof),0,atol=3e-15)

    # Independent slow scatter reference on a two-triangle patch.
    reference=np.zeros(K.shape)
    for e,ids in enumerate(mesh.elements):
        local=construct_k_matrix(mesh.p1_points[mesh.p1_elements[e]],quad)
        for i,a in enumerate(ids):
            for j,b in enumerate(ids):
                reference[a,b]+=local[i,j]
    np.testing.assert_allclose(K.toarray(),reference,atol=2e-15)

    # Exact constant load, including extrapolated flux and signed spans.
    ve=np.array([[1.5,0.],[3.5,.2],[1.8,1.]])
    J=np.column_stack((ve[1]-ve[0],ve[2]-ve[0]))
    physical_nodes=P3_REFERENCE_NODES@J.T+ve[0]
    grid=np.linspace(0,1,9)
    s0=3.5
    for span in (-.7,.7):
        psi_e=span*(physical_nodes[:,0]-3.)
        pp=np.zeros_like(grid)
        ff=span*s0*(2.+grid)
        fe=construct_rhs(ve,psi_e,quad,pp,ff,span)
        np.testing.assert_allclose(fe.sum(),s0*np.linalg.det(J)/2,rtol=2e-15,atol=2e-15)

    # Pressure term and piecewise linear interpolation/extrapolation.
    psi_e=np.linspace(-1.3,.7,10)
    span=-.8
    pp=np.array([-2e4,-3e4,-1e4,-4e4])
    ff=np.array([.2,-.1,.4,.3])
    u=np.linspace(0,1,4)
    query=(quad.N@psi_e+span)/span
    index=np.clip(np.searchsorted(u,query,side="right")-1,0,len(u)-2)
    def linear_reference(values):
        slope=(values[index+1]-values[index])/(u[index+1]-u[index])
        return values[index]+slope*(query-u[index])
    rq=ve[0,0]+J[0,0]*quad.xi+J[0,1]*quad.eta
    source=(MU0*rq*linear_reference(pp)+linear_reference(ff)/rq)/span
    expected=quad.N.T@(quad.w*np.linalg.det(J)*source)
    actual=construct_rhs(ve,psi_e,quad,pp,ff,span)
    np.testing.assert_allclose(actual,expected,rtol=5e-14,atol=2e-16)
    np.testing.assert_allclose(
        construct_rhs(ve,-psi_e,quad,pp,ff,-span),-actual,rtol=5e-14,atol=2e-16)

    # Constant profiles are independent of the supplied flux at fixed span.
    a=construct_rhs(ve,np.arange(10.),quad,np.ones(4),np.ones(4),-.5)
    b=construct_rhs(ve,np.zeros(10),quad,np.ones(4),np.ones(4),-.5)
    np.testing.assert_allclose(a,b,rtol=0,atol=0)

    psi=np.linspace(-.8,.1,mesh.n_dof)
    with patch.object(fem,"interp1d",wraps=fem.interp1d) as interpolation:
        global_rhs=assemble_rhs(mesh,psi,quad,pp,ff,span)
        assert interpolation.call_count==1
    reference_rhs=np.zeros(mesh.n_dof)
    for e,ids in enumerate(mesh.elements):
        local=construct_rhs(mesh.p1_points[mesh.p1_elements[e]],psi[ids],quad,pp,ff,span)
        np.add.at(reference_rhs,ids,local)
    np.testing.assert_allclose(global_rhs,reference_rhs,rtol=2e-15,atol=2e-16)
    assert not np.allclose(global_rhs,assemble_rhs(mesh,psi*.3,quad,pp,ff,span))

    # Nonzero boundary reduction and exact vacuum field.
    exact=mesh.points[:,0]**2
    zero=np.zeros(mesh.n_dof)
    for boundary in (exact,exact[mesh.boundary_nodes],lambda R,Z:R**2):
        KII,rhsI,bc=apply_bc(mesh,K,zero,boundary)
        np.linalg.cholesky(KII.toarray())
        result=bc.restore(spsolve(KII,rhsI))
        np.testing.assert_allclose(result,exact,rtol=2e-14,atol=2e-14)
        np.testing.assert_allclose(bc.reduce_rhs(global_rhs),global_rhs[bc.free_nodes]-bc.boundary_correction)
    for matrix in (K.tocoo(),K.tolil(),K.todok()):
        KII,rhsI,bc=apply_bc(mesh,matrix,zero,.37)
        np.testing.assert_allclose(bc.restore(spsolve(KII,rhsI)),.37,rtol=2e-14,atol=2e-14)
    KII,rhsI,bc=apply_bc(mesh,K,zero)
    solve=factorized(KII.tocsc())
    np.testing.assert_allclose(solve(bc.reduce_rhs(zero)),0)

    # A second vacuum field has a nonpolynomial integrand: verify quadrature convergence.
    exact_mixed=mesh.points[:,0]**2+mesh.points[:,1]
    quadrature_errors=[]
    for degree in (10,20):
        A=assemble_k_matrix(mesh,precompute_quadrature(degree))
        Aii,fi,partition=apply_bc(mesh,A,zero,exact_mixed)
        quadrature_errors.append(np.max(abs(partition.restore(spsolve(Aii,fi))-exact_mixed)))
    assert quadrature_errors[1]<1e-13 and quadrature_errors[1]<quadrature_errors[0]/100

    # Actual concave-domain solve with the exact vacuum field R^2.
    notch=np.array([[2.,0.],[5.,0.],[5.,1.],[3.,1.],[3.,2.],[5.,2.],[5.,3.],[2.,3.]])
    concave=generate_p3_nodes(build_p1_mesh(notch,nr=13,nz=15))
    matrix=assemble_k_matrix(concave,quad)
    exact=concave.points[:,0]**2
    reduced,rhsI,bc=apply_bc(concave,matrix,np.zeros(concave.n_dof),exact)
    error=np.max(abs(bc.restore(spsolve(reduced,rhsI))-exact))
    assert error<1e-10

    # Invalid inputs.
    def fails(fn,*args,**kwargs):
        try: fn(*args,**kwargs)
        except (ValueError,TypeError): return
        raise AssertionError("Invalid input accepted.")
    for bad in (0.,np.nan,np.inf,True,[1.],1j):
        fails(assemble_rhs,mesh,psi,quad,pp,ff,bad)
    fails(assemble_rhs,mesh,psi[:-1],quad,pp,ff,span)
    fails(assemble_rhs,mesh,psi+1j,quad,pp,ff,span)
    fails(assemble_rhs,mesh,psi,quad,pp,ff[:-1],span)
    fails(construct_rhs,ve,psi_e[:-1],quad,pp,ff,span)
    fails(apply_bc,mesh,K,zero,[1.,2.])
    fails(apply_bc,mesh,K,zero,np.nan)
    fails(apply_bc,mesh,K.astype(complex)+1j*K,zero)
    fails(apply_bc,mesh,K.toarray(),zero)
    fails(apply_bc,replace(mesh,boundary_nodes=np.array([0,0])),K,zero)
    fails(apply_bc,replace(mesh,boundary_nodes=np.arange(mesh.n_dof)),K,zero)
    fails(assemble_k_matrix,replace(mesh,elements=np.zeros_like(mesh.elements)),quad)
    print("Assembly, normalized-profile RHS, extrapolation, signed-span, and Dirichlet checks passed.")
    print("Concave-domain R^2 solve:",concave.n_dof,"DOFs, max error",error)


if __name__ == "__main__":
    test_assembly()
