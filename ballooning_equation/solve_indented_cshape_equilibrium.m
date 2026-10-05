function [equilibrium, diagnostics, predictor] = solve_indented_cshape_equilibrium( ...
    config, shape, p0, prepared, warm)
%SOLVE_INDENTED_CSHAPE_EQUILIBRIUM Recover failed pressure steps by continuation.
%   Only equilibria satisfying both GS and prescribed-q tolerances can be
%   predictors. Auxiliary pressures do not become stability sampling points.
    if nargin<4 || isempty(prepared)
        input=build_indented_cshape_scan_input(config,shape,p0);
        prepared=prepare_GS_problem(input);
    end
    if nargin<5, warm=[]; end
    if ~isfield(prepared.qEvaluation,'contourCache')
        prepared.qEvaluation.contourCache=build_q_contour_cache( ...
            prepared.mesh,prepared.qEvaluation.contourNR);
    end
    controls=struct('Enabled',false, ...
        'BootstrapPressureAmplitude',.002*config.BReference^2/(2*(4*pi*1.e-7)), ...
        'MinimumAbsoluteStep',10,'MinimumRelativeStep',1.e-3, ...
        'StepGrowth',1.5,'MaxAttempts',24,'MaximumQIterationRestarts',1);
    if isfield(config,'Continuation'), controls=config.Continuation; end
    if ~isfield(controls,'MaximumQIterationRestarts'), controls.MaximumQIterationRestarts=1; end
    validateattributes(controls.Enabled,{'logical'},{'scalar'});
    validateattributes(controls.BootstrapPressureAmplitude,{'double'},{'scalar','finite','positive'});
    validateattributes(controls.MinimumAbsoluteStep,{'double'},{'scalar','finite','nonnegative'});
    validateattributes(controls.MinimumRelativeStep,{'double'},{'scalar','finite','positive'});
    validateattributes(controls.StepGrowth,{'double'},{'scalar','finite','>',1});
    validateattributes(controls.MaxAttempts,{'double'},{'scalar','integer','>=',2});
    validateattributes(controls.MaximumQIterationRestarts,{'double'},{'scalar','integer','nonnegative'});
    if ~isempty(warm) && (~warm.solverResult.converged || ~warm.constraintResult.converged)
        error('cshape:Continuation:InvalidPredictor','The predictor must satisfy GS and q convergence.');
    end
    predictor=warm;
    diagnostics=struct('attempts',struct([]),'usedContinuation',false, ...
        'lastConvergedPressure',NaN,'terminationReason','direct_failure','lastFailure',[]);
    [equilibrium,record]=attempt(p0,warm,'target');
    diagnostics.attempts=record;
    if record.success
        predictor=equilibrium;
        diagnostics.lastConvergedPressure=p0;
        diagnostics.terminationReason='converged'; return
    end
    diagnostics.lastFailure=record.failure;
    if ~controls.Enabled, return; end
    diagnostics.usedContinuation=true;
    if isempty(warm)
        bootstrap=min(controls.BootstrapPressureAmplitude,p0/10);
        if bootstrap<=0
            diagnostics.terminationReason='no_bootstrap_pressure'; return
        end
        [predictor,record]=attempt(bootstrap,[],'bootstrap');
        diagnostics.attempts(end+1)=record;
        if ~record.success
            diagnostics.lastFailure=record.failure;
            diagnostics.terminationReason='bootstrap_failure'; return
        end
    end
    current=pressure_of(predictor);
    diagnostics.lastConvergedPressure=current;
    if abs(p0-current)<=100*eps(max(1,abs(p0)))
        diagnostics.terminationReason='same_pressure_predictor_failed'; return
    end
    step=abs(p0-current)/2;
    direction=sign(p0-current);
    minimumStep=max(controls.MinimumAbsoluteStep,controls.MinimumRelativeStep*abs(p0));
    while numel(diagnostics.attempts)<controls.MaxAttempts
        remaining=abs(p0-current);
        if remaining<=100*eps(max(1,abs(p0)))
            diagnostics.terminationReason='unresolved_target'; return
        end
        if step<minimumStep && step<remaining
            diagnostics.terminationReason='minimum_pressure_step'; return
        end
        trialPressure=current+direction*min(step,remaining);
        if step>=remaining, trialPressure=p0; end
        [trial,record]=attempt(trialPressure,predictor,'pressure_step');
        diagnostics.attempts(end+1)=record;
        if record.success
            predictor=trial; current=trialPressure;
            diagnostics.lastConvergedPressure=current;
            if trialPressure==p0
                equilibrium=trial;
                diagnostics.terminationReason='converged'; return
            end
            step=step*controls.StepGrowth;
        else
            diagnostics.lastFailure=record.failure;
            % The target can clip a grown step. Halve the pressure gap
            % actually attempted, rather than repeating the identical
            % clipped target with the same predictor several times.
            step=min(step,remaining)/2;
        end
    end
    diagnostics.terminationReason='maximum_attempts';

    function [eq,record]=attempt(pressure,initial,stage)
        timer=tic;
        record=struct('pressureAmplitude',pressure,'stage',stage,'success',false, ...
            'innerResidual',NaN,'qMaximumRelativeError',NaN, ...
            'qIterationRestarts',0,'wallTime',NaN,'failure',[]);
        try
            input=build_indented_cshape_scan_input(config,shape,pressure,initial);
            problem=prepared;
            if ~isempty(initial) && isequal(initial.mesh.P1points,problem.mesh.P1points) && ...
                    isequal(initial.mesh.elements,problem.mesh.elements)
                problem.initial.psi=initial.psi;
            end
            eq=solve_GS(input,problem);
            % Continue the same-pressure q optimization after an iteration
            % budget is exhausted. This unfinished q iterate is never a
            % predictor at another pressure or a stability sampling point.
            while eq.solverResult.converged && ~eq.constraintResult.converged && ...
                    strcmp(eq.constraintResult.terminationReason,'maximum_iterations') && ...
                    record.qIterationRestarts<controls.MaximumQIterationRestarts
                if config.Verbose
                    fprintf('  q iteration restart at p0 %.6g, q error %.3e\n', ...
                        pressure,eq.constraintResult.finalMaxAbsRelativeQ);
                end
                input=build_indented_cshape_scan_input(config,shape,pressure,eq);
                problem.initial.psi=eq.psi;
                eq=solve_GS(input,problem);
                record.qIterationRestarts=record.qIterationRestarts+1;
            end
            record.innerResidual=eq.solverResult.finalRelativeResidual;
            record.qMaximumRelativeError=eq.constraintResult.finalMaxAbsRelativeQ;
            if ~eq.solverResult.converged || ~eq.constraintResult.converged
                error('cshape:Point:Nonconvergence', ...
                    'GS or q constraint did not converge (%s); residual %.3e, q error %.3e.', ...
                    eq.constraintResult.terminationReason,record.innerResidual,record.qMaximumRelativeError);
            end
            eq.scanPressureAmplitude=pressure;
            record.success=true;
        catch exception
            if strcmp(exception.identifier,'MATLAB:OperationTerminated'), rethrow(exception); end
            eq=[];
            record.failure=struct('identifier',exception.identifier,'message',exception.message);
        end
        record.wallTime=toc(timer);
        if config.Verbose && (~strcmp(stage,'target') || ~record.success)
            fprintf('  continuation p0 %.6g: converged %d, residual %.3e, q error %.3e (%.1fs)\n', ...
                pressure,record.success,record.innerResidual,record.qMaximumRelativeError,record.wallTime);
        end
    end
end

function pressure=pressure_of(eq)
    if isfield(eq,'scanPressureAmplitude')
        pressure=eq.scanPressureAmplitude;
    else
        u=eq.profile.psiN(:);
        basis=-4*u.*(1-u.^2);
        pressure=(basis.'*eq.profile.dp_dpsiN(:))/(basis.'*basis);
    end
end
