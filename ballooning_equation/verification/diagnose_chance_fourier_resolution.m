function result=diagnose_chance_fourier_resolution(equilibrium,outputFile)
%DIAGNOSE_CHANCE_FOURIER_RESOLUTION Separate spectral truncation from GS noise.
    setup_indented_cshape_scan;
    [eq,fun]=build_ballooning_adapter(equilibrium,NR=1000);
    labels=[.85 .9 .95]; modes=[16 24 32 40 48 64 96];
    result=struct('psiN',labels,'fourierModes',modes,'lambda',nan(3,numel(modes)), ...
        'closure',nan(3,numel(modes)),'curvatureMaximum',nan(3,numel(modes)), ...
        'fluxRms',nan(3,numel(modes)),'qCorrection',nan(3,numel(modes)), ...
        'cellsPerTurn',384,'nSurfacePoints',8192,'failures',{cell(3,numel(modes))});
    for k=1:numel(labels)
        for j=1:numel(modes)
            try
                state=prepare_salpha_equilibrium(eq,fun,labels(k),Model='Direct', ...
                    NSurfacePoints=8192,FourierModes=modes(j));
                coef=evaluate_gcf(state.thetaDiag,0,theta_bnd=5*pi, ...
                    a_N=equilibrium.reference.aMinor,B_N=abs(equilibrium.reference.B0));
                theta=linspace(-5*pi,5*pi,5*384+1).';
                [G,C,F,t]=construct_matrix(theta,interp1(coef.theta,coef.g,theta,'pchip'), ...
                    interp1(coef.theta,coef.c,theta,'pchip'),interp1(coef.theta,coef.f,theta,'pchip'));
                result.lambda(k,j)=calculate_eigenvalue(t,G,C,F,'dirichlet',NumEigenvalues=1);
                d=state.diagnostics.directGeometry;
                result.closure(k,j)=state.thetaDiag.input_closure_defect;
                result.curvatureMaximum(k,j)=max(abs(state.thetaDiag.curvature_normal));
                result.fluxRms(k,j)=d.psiNRms; result.qCorrection(k,j)=d.qCorrection;
                fprintf('u %.3f modes %2d: lambda %+.8e, closure %.4g, max curvature %.5g\n', ...
                    labels(k),modes(j),result.lambda(k,j),result.closure(k,j),result.curvatureMaximum(k,j));
            catch problem
                result.failures{k,j}=struct('identifier',problem.identifier,'message',problem.message);
                fprintf('u %.3f modes %d: %s\n',labels(k),modes(j),problem.message);
            end
            save(outputFile,'result','-v7.3');
        end
    end
end
