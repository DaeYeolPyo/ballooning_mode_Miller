function result = evaluate_equilibrium_ballooning(equilibrium, opts)
%EVALUATE_EQUILIBRIUM_BALLOONING Actual Direct geometry, sampled u and phase.
%   A positive eigenvalue proves instability even if another solve fails.
%   lambdaMax is NaN for incomplete scans; lambdaMaxObserved remains a
%   lower bound. An incomplete scan can never be declared stable.
    arguments
        equilibrium (1,1) struct
        opts.RadialLabels double = [.05:.1:.85,.90,.95]
        opts.Theta0 double = 2*pi*((0:7)-4)/8
        opts.ThetaBoundary (1,1) double {mustBeFinite,mustBePositive} = 5*pi
        opts.CellsPerTurn (1,1) double {mustBeFinite,mustBeInteger,mustBeGreaterThanOrEqual(opts.CellsPerTurn,8)} = 24
        opts.NSurfacePoints (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.NSurfacePoints,32)} = 1024
        opts.FourierModes (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.FourierModes,2)} = 48
        opts.ProjectionIterations (1,1) double {mustBeInteger,mustBeNonnegative} = 2
        opts.AdapterNR (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(opts.AdapterNR,16)} = 257
        opts.EnforceDirectQ (1,1) logical = true
        opts.EnforcePeriodicity (1,1) logical = true
        opts.LambdaTolerance (1,1) double {mustBeFinite,mustBeNonnegative} = 1.e-6
        opts.MaximumFluxRms (1,1) double {mustBeFinite,mustBePositive} = 2.e-3
        opts.MaximumFluxError (1,1) double {mustBeFinite,mustBePositive} = 1.e-2
        opts.MinimumNormalAlignment (1,1) double {mustBeFinite} = .95
        opts.MaximumQCorrection (1,1) double {mustBeFinite,mustBePositive} = .05
        opts.MaximumInputClosureDefect (1,1) double {mustBeFinite,mustBePositive} = .05
        opts.StoreWorstEigenfunction (1,1) logical = true
        opts.Verbose (1,1) logical = false
    end
    setup_indented_cshape_scan;
    validateattributes(opts.RadialLabels,{'double'},{'real','vector','nonempty','finite','>',0,'<',1});
    validateattributes(opts.Theta0,{'double'},{'real','vector','nonempty','finite'});
    u=sort(opts.RadialLabels(:));
    phase=opts.Theta0(:).';
    if numel(unique(u))~=numel(u) || ...
            numel(unique(round(mod(phase,2*pi)*1.e12)))~=numel(phase)
        error('cshape:Ballooning:RepeatedGrid','Require distinct surfaces and phases modulo 2*pi.');
    end
    if opts.FourierModes>floor((opts.NSurfacePoints-1)/2) || ...
            opts.MinimumNormalAlignment<0 || opts.MinimumNormalAlignment>1
        error('cshape:Ballooning:Options','Invalid Fourier resolution or normal alignment threshold.');
    end
    if isfield(equilibrium,'solverResult') && ~equilibrium.solverResult.converged
        error('cshape:Ballooning:GSConvergence','Refuse an unconverged GS equilibrium.');
    end
    if isfield(equilibrium,'constraintResult') && ~equilibrium.constraintResult.converged
        error('cshape:Ballooning:QConvergence','Refuse an unconverged q constraint.');
    end
    [eq,eqfunc,adapterDiagnostics]=build_ballooning_adapter(equilibrium,NR=opts.AdapterNR);
    aN=equilibrium.reference.aMinor;
    BN=abs(equilibrium.reference.B0);
    ncell=2*ceil(opts.ThetaBoundary*opts.CellsPerTurn/(2*pi));
    theta=linspace(-opts.ThetaBoundary,opts.ThetaBoundary,ncell+1).';
    values=nan(numel(u),numel(phase));
    errors=cell(numel(u),numel(phase));
    surfaceDiagnostics=cell(numel(u),1);
    worst=struct('psiN',NaN,'theta0',NaN,'lambda',-inf,'theta',[],'X',[]);
    for j=1:numel(u)
        try
            state=prepare_salpha_equilibrium(eq,eqfunc,u(j),Model="Direct", ...
                NSurfacePoints=opts.NSurfacePoints,FourierModes=opts.FourierModes, ...
                ProjectionIterations=opts.ProjectionIterations, ...
                EnforceDirectQ=opts.EnforceDirectQ,EnforcePeriodicity=opts.EnforcePeriodicity);
            d=state.diagnostics.directGeometry;
            surfaceDiagnostics{j}=struct('geometry',d,'thetaPsi',struct( ...
                'inputClosure',state.thetaDiag.closure_input, ...
                'usedClosure',state.thetaDiag.closure_used, ...
                'inputClosureDefect',state.thetaDiag.input_closure_defect, ...
                'FFprimeInput',state.thetaDiag.FFprime_input, ...
                'FFprimeUsed',state.thetaDiag.FFprime_used));
            quality=[d.psiNRms,d.psiNMaximum,d.minimumNormalAlignment, ...
                d.qCorrection,state.thetaDiag.input_closure_defect];
            if any(~isfinite(quality)) || ...
                    d.psiNRms>opts.MaximumFluxRms || d.psiNMaximum>opts.MaximumFluxError || ...
                    d.minimumNormalAlignment<opts.MinimumNormalAlignment || ...
                    abs(d.qCorrection-1)>opts.MaximumQCorrection || ...
                    state.thetaDiag.input_closure_defect>opts.MaximumInputClosureDefect
                error('cshape:Ballooning:GeometryQuality', ...
                    ['psiN %.4g: surface quality failed; psi RMS %.3e (limit %.3e), ', ...
                     'max %.3e (%.3e), normal %.6f (min %.6f), ', ...
                     'q correction %.6f (deviation limit %.3e), closure %.3e (%.3e).'], ...
                    u(j),d.psiNRms,opts.MaximumFluxRms,d.psiNMaximum,opts.MaximumFluxError, ...
                    d.minimumNormalAlignment,opts.MinimumNormalAlignment, ...
                    d.qCorrection,opts.MaximumQCorrection, ...
                    state.thetaDiag.input_closure_defect,opts.MaximumInputClosureDefect);
            end
        catch exception
            for k=1:numel(phase), errors{j,k}=failure(exception); end
            continue
        end
        for k=1:numel(phase)
            try
                coeff=evaluate_gcf(state.thetaDiag,phase(k), ...
                    theta_bnd=opts.ThetaBoundary,a_N=aN,B_N=BN);
                g=interp1(coeff.theta,coeff.g,theta,'pchip');
                c=interp1(coeff.theta,coeff.c,theta,'pchip');
                f=interp1(coeff.theta,coeff.f,theta,'pchip');
                if any(g<=0) || any(f<=0) || any(~isfinite([g;c;f]))
                    error('cshape:Ballooning:Coefficients','Invalid g, c or f.');
                end
                [G,C,F,tdof]=construct_matrix(theta,g,c,f);
                [lambda,X]=calculate_eigenvalue(tdof,G,C,F,'dirichlet',NumEigenvalues=1);
                if isempty(lambda) || ~isfinite(lambda(1))
                    error('cshape:Ballooning:Eigenvalue','No finite largest eigenvalue.');
                end
                values(j,k)=lambda(1);
                if lambda(1)>worst.lambda
                    worst.psiN=u(j); worst.theta0=phase(k); worst.lambda=lambda(1);
                    if opts.StoreWorstEigenfunction
                        worst.theta=tdof; worst.X=X(:,1)/max(abs(X(:,1)));
                    end
                end
            catch exception
                errors{j,k}=failure(exception);
            end
        end
        if opts.Verbose
            fprintf('  psiN %.4f: %d/%d phase solves\n',u(j),sum(isfinite(values(j,:))),numel(phase));
        end
    end
    complete=all(isfinite(values),'all');
    observed=values(isfinite(values));
    if isempty(observed), observedMaximum=NaN; else, observedMaximum=max(observed); end
    status="unknown"; lambdaMax=NaN;
    if complete, lambdaMax=observedMaximum; end
    if observedMaximum>opts.LambdaTolerance
        status="unstable";
    elseif complete && lambdaMax < -opts.LambdaTolerance
        status="stable";
    elseif complete
        status="marginal";
    end
    if isempty(observed), worst.lambda=NaN; end
    result=struct('status',status,'complete',complete,'lambdaMax',lambdaMax, ...
        'lambdaMaxObserved',observedMaximum,'lambda',values,'psiN',u,'theta0',phase, ...
        'surfaceLambdaMax',max(values,[],2,'omitnan'),'worst',worst, ...
        'failures',{errors},'numberFailed',sum(~isfinite(values),'all'), ...
        'surfaceDiagnostics',{surfaceDiagnostics},'adapterDiagnostics',adapterDiagnostics, ...
        'normalization',struct('aN',aN,'BN',BN),'options',opts);
end

function item=failure(exception)
    item=struct('identifier',exception.identifier,'message',exception.message);
end
