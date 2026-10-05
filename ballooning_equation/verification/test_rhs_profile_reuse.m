function report = test_rhs_profile_reuse()
%TEST_RHS_PROFILE_REUSE Compare cached and original element assembly.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(root); setup_indented_cshape_scan;
    saved=load(fullfile(root,'output','indented_cshape_dense','column_03.mat'),'scan');
    eq=saved.scan.equilibria{find(saved.scan.done,1,'last')};
    P3=eq.mesh; psi=eq.psi;
    p=eq.profile.dp_dpsiN; ff=eq.profile.FdF_dpsiN;
    dpsi=eq.global.sibry-eq.global.simag;
    maximumError=0;
    for degree=[7,10]
        quad=precompute_quadrature(int8(degree));
        for orientation=1:2
            if orientation==2, p=p.'; ff=ff.'; end
            for fluxSign=[1,-1]
                original=reference_rhs(P3,fluxSign*psi,quad,p,ff,fluxSign*dpsi);
                reused=assemble_RHS(P3,fluxSign*psi,quad,p,ff,fluxSign*dpsi);
                err=norm(reused-original,inf)/max(1,norm(original,inf));
                maximumError=max(maximumError,err);
                assert(err<1.e-14,'Cached RHS differs from original element assembly.');
            end
        end
    end
    quad=precompute_quadrature(int8(7));
    K=assemble_K_matrix(P3,quad,false);
    symmetryError=norm(K-K.',inf)/max(1,norm(K,inf));
    assert(symmetryError<1.e-13);
    reference_rhs(P3,psi,quad,p,ff,dpsi); assemble_RHS(P3,psi,quad,p,ff,dpsi);
    repetitions=5; timer=tic;
    for j=1:repetitions, reference_rhs(P3,psi,quad,p,ff,dpsi); end
    originalTime=toc(timer)/repetitions; timer=tic;
    for j=1:repetitions, assemble_RHS(P3,psi,quad,p,ff,dpsi); end
    reusedTime=toc(timer)/repetitions;
    report=struct('maximumRelativeRHSError',maximumError,'stiffnessSymmetryError',symmetryError, ...
        'originalSecondsPerAssembly',originalTime,'reusedSecondsPerAssembly',reusedTime, ...
        'assemblySpeedup',originalTime/reusedTime,'numberElements',size(P3.elements,1));
    fprintf('RHS agreement %.3g; stiffness symmetry %.3g; original %.4fs reused %.4fs (%.2fx)\n', ...
        maximumError,symmetryError,originalTime,reusedTime,report.assemblySpeedup);
end

function f=reference_rhs(P3,psi,quad,p,ff,dpsi)
    f=zeros(size(P3.points,1),1);
    for e=1:size(P3.elements,1)
        ids=P3.elements(e,:);
        V=P3.P1points(P3.P1elements(e,:),:);
        fe=construct_RHS(V,psi(ids),quad,p,ff,dpsi);
        f(ids)=f(ids)+fe;
    end
end
