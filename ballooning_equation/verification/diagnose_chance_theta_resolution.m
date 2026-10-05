function result=diagnose_chance_theta_resolution(equilibrium,outputFile)
%DIAGNOSE_CHANCE_THETA_RESOLUTION Hold native equilibrium/geometry fixed.
%   Bean shaping compresses the bad-curvature region in PEST theta. Resolve
%   that coefficient structure rather than trusting a coarse FE eigenvalue.
    setup_indented_cshape_scan;
    [eq,fun]=build_ballooning_adapter(equilibrium,NR=1000);
    labels=[.9 .95]; cells=[24 48 96 192 384 768 1536];
    values=nan(numel(labels),numel(cells)); quality=cell(numel(labels),1);
    result=struct('psiN',labels,'cellsPerTurn',cells,'lambda',values, ...
        'quality',{quality},'thetaBoundary',5*pi,'nSurfacePoints',8192,'fourierModes',64);
    for k=1:numel(labels)
        state=prepare_salpha_equilibrium(eq,fun,labels(k),Model='Direct', ...
            NSurfacePoints=8192,FourierModes=64);
        quality{k}=struct('geometry',state.diagnostics.directGeometry, ...
            'closure',state.thetaDiag.input_closure_defect);
        coef=evaluate_gcf(state.thetaDiag,0,theta_bnd=5*pi, ...
            a_N=equilibrium.reference.aMinor,B_N=abs(equilibrium.reference.B0));
        for j=1:numel(cells)
            n=2*ceil(5*pi*cells(j)/(2*pi)); theta=linspace(-5*pi,5*pi,n+1).';
            [G,C,F,t]=construct_matrix(theta,interp1(coef.theta,coef.g,theta,'pchip'), ...
                interp1(coef.theta,coef.c,theta,'pchip'),interp1(coef.theta,coef.f,theta,'pchip'));
            values(k,j)=calculate_eigenvalue(t,G,C,F,'dirichlet',NumEigenvalues=1);
            fprintf('u %.3f, theta cells/turn %4d: lambda %+.9e\n',labels(k),cells(j),values(k,j));
            result.lambda=values; result.quality=quality;
            save(outputFile,'result','-v7.3');
        end
    end
end
