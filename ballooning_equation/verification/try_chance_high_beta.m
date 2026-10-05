function result = try_chance_high_beta(checkpoint, outputFile, opts)
%TRY_CHANCE_HIGH_BETA Diagnose a physical pressure step; never classify bad q.
    arguments
        checkpoint (1,1) string
        outputFile (1,1) string
        opts.PressureAmplitude (1,1) double = 1.04615e6
        opts.GridSize (1,2) double = [48 64]
        opts.ContourNR (1,1) double = 800
        opts.PressureCompensation (1,1) logical = true
        opts.RelativeFD (1,1) double = .02
        opts.NBoundary (1,1) double = 128
    end
    setup_indented_cshape_scan;
    saved=load(checkpoint,'report'); source=saved.report;
    candidates=find(~cellfun(@isempty,source.equilibria));
    warm=source.equilibria{candidates(end)}; shape=source.shapes{1};
    config=make_chance_prl_figure3_config(Indentation=shape.iIn,GridSize=opts.GridSize);
    config.Grid.nBoundary=opts.NBoundary;
    shape=chance_prl_bean_boundary(shape.iIn,NBoundary=opts.NBoundary);
    config.QEvaluation.contourNR=opts.ContourNR;
    config.Ballooning.AdapterNR=opts.ContourNR;
    config.PQConstraint.FFprime.finiteDifferenceRelativeStep=opts.RelativeFD;
    config.PQConstraint.FFprime.finiteDifferenceAbsoluteStep=.05*opts.RelativeFD/.02;
    input=build_indented_cshape_scan_input(config,shape,opts.PressureAmplitude,warm);
    if opts.PressureCompensation
        u=warm.profile.psiN;
        newDp=-4*opts.PressureAmplitude*u.*(1-u.^2);
        deltaDp=newDp-warm.profile.dp_dpsiN;
        % Cylindrical source balance is only an initial prediction. The full
        % GS/q solve below must still meet its original residual tolerances.
        input.profile.initial.FFprime=warm.profile.FdF_dpsiN ...
            -(4*pi*1.e-7)*shape.Rgeo^2*deltaDp;
    end
    timer=tic;
    prepared=prepare_GS_problem(input);
    prepared.qEvaluation.contourCache=build_q_contour_cache(prepared.mesh,opts.ContourNR);
    oldMesh=triangulation(warm.mesh.P1elements,warm.mesh.P1points);
    projected=eval_P3_sol(oldMesh,warm.mesh,warm.psi,prepared.mesh.points);
    if all(isfinite(projected)), prepared.initial.psi=projected; end
    result=struct('config',config,'shape',shape,'options',opts,'equilibrium',[], ...
        'beta',[],'ballooning',[],'failure',[],'wallTime',NaN);
    try
        equilibrium=solve_GS(input,prepared);
        result.equilibrium=equilibrium;
        fprintf('p0 %.7g: %s, GS %.3e, q max %.3e\n',opts.PressureAmplitude, ...
            equilibrium.constraintResult.terminationReason, ...
            equilibrium.solverResult.finalRelativeResidual, ...
            equilibrium.constraintResult.finalMaxAbsRelativeQ);
        if equilibrium.solverResult.converged && equilibrium.constraintResult.converged
            equilibrium.scanPressureAmplitude=opts.PressureAmplitude;
            result.equilibrium=equilibrium;
            result.beta=compute_volume_beta(equilibrium);
            args=namedargs2cell(config.Ballooning);
            result.ballooning=evaluate_equilibrium_ballooning(equilibrium,args{:});
            fprintf('beta %.7g%%, %s, lambda %+.6e, complete %d\n', ...
                result.beta.betaPercent,result.ballooning.status, ...
                result.ballooning.lambdaMax,result.ballooning.complete);
        end
    catch problem
        result.failure=struct('identifier',problem.identifier,'message',problem.message);
        fprintf('%s: %s\n',problem.identifier,problem.message);
    end
    result.wallTime=toc(timer);
    save(outputFile,'result','-v7.3');
end
