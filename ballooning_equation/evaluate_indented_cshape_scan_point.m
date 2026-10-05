function [point, equilibrium, predictor] = evaluate_indented_cshape_scan_point(config, shape, p0, prepared, warm)
%EVALUATE_INDENTED_CSHAPE_SCAN_POINT One self-consistent equilibrium.
    if nargin<4, prepared=[]; end
    if nargin<5, warm=[]; end
    timer=tic;
    equilibrium=[];
    predictor=warm;
    point=struct('aspectRatio',shape.aspectRatio,'indent',shape.indent, ...
        'iOut',shape.iOut,'pressureAmplitude',p0,'betaPercent',NaN, ...
        'status',"unknown",'complete',false,'equilibriumConverged',false, ...
        'lambdaMax',NaN,'lambdaMaxObserved',NaN,'betaDiagnostics',[], ...
        'ballooning',[],'solverDiagnostics',[],'equilibriumAttempts',[], ...
        'failure',[],'wallTime',NaN);
    stage='equilibrium';
    try
        [equilibrium,attempts,predictor]=solve_indented_cshape_equilibrium(config,shape,p0,prepared,warm);
        point.equilibriumAttempts=attempts;
        if isempty(equilibrium)
            error('cshape:Point:EquilibriumFailure','%s. Last failure: %s', ...
                attempts.terminationReason,attempts.lastFailure.message);
        end
        point.equilibriumConverged=true;
        q=equilibrium.qResult.q(:);
        qu=equilibrium.qResult.psiN(:);
        qt=polyval(fliplr(config.QCoefficients(:).'),qu);
        point.solverDiagnostics=struct( ...
            'innerResidual',equilibrium.solverResult.finalRelativeResidual, ...
            'innerIterations',equilibrium.solverResult.iterations, ...
            'qTerminationReason',equilibrium.constraintResult.terminationReason, ...
            'qOuterIterations',equilibrium.constraintResult.outerIterations, ...
            'qConstraintMaximumRelativeError',equilibrium.constraintResult.finalMaxAbsRelativeQ, ...
            'qConstraintRmsLogError',equilibrium.constraintResult.finalRmsLogQ, ...
            'qMaximumRelativeError',max(abs((q-qt)./qt)), ...
            'Ip',equilibrium.global.Ip,'Raxis',equilibrium.global.Raxis, ...
            'nativeFluxSpan',equilibrium.global.sibry-equilibrium.global.simag);
        stage='beta';
        beta=compute_volume_beta(equilibrium,QuadratureDegree=config.BetaQuadratureDegree);
        point.betaDiagnostics=beta;
        point.betaPercent=beta.betaPercent;
        stage='ballooning';
        argumentsCell=namedargs2cell(config.Ballooning);
        bal=evaluate_equilibrium_ballooning(equilibrium,argumentsCell{:});
        point.ballooning=bal;
        point.lambdaMax=bal.lambdaMax;
        point.lambdaMaxObserved=bal.lambdaMaxObserved;
        point.complete=bal.complete;
        point.status=bal.status;
    catch exception
        if strcmp(exception.identifier,'MATLAB:OperationTerminated'), rethrow(exception); end
        point.failure=struct('stage',stage,'identifier',exception.identifier,'message',exception.message);
        if ~point.equilibriumConverged, equilibrium=[]; end
    end
    point.wallTime=toc(timer);
end
