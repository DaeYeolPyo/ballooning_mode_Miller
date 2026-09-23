function [solution, result] = solve_pq_constraint( ...
    problem, pprime, qTarget, FFprimeModel0, psi0, options)
%SOLVE_PQ_CONSTRAINT Solve prescribed-p' prescribed-q equilibria.
%
% The unknowns are low-order coefficients of FF'(psiN). Every residual
% and every accepted step is evaluated by a complete inner GS solve. The
% old frozen-geometry log(G) update is deliberately not used.

    narginchk(5, 6);

    if nargin < 6 || isempty(options)
        options = struct();
    end

    assert_required_fields(problem, ...
        {'mesh','profile','initial','qEvaluation'}, 'problem');
    assert_required_fields(problem.profile, {'psiN'}, 'problem.profile');

    if isempty(psi0)
        psi0 = problem.initial.psi;
    end

    profilePsiN = problem.profile.psiN(:);
    pprime = pprime(:);
    psi0 = psi0(:);

    if numel(pprime) ~= numel(profilePsiN) || any(~isfinite(pprime))
        error('solve_pq_constraint:InvalidPprime', ...
            'pprime must be finite and sampled on problem.profile.psiN.');
    end

    assert_required_fields(FFprimeModel0, ...
        {'coefficients','nBasis','Fboundary','signF'}, ...
        'FFprimeModel0');

    outerOptions = get_substructure(options, 'outer');
    predictorOptions = get_substructure(options, 'FFprime');

    if isfield(options,'inner') && isstruct(options.inner)
        innerOptions = options.inner;
    elseif isfield(problem,'options') && ...
            isfield(problem.options,'Picard')
        innerOptions = problem.options.Picard;
    else
        innerOptions = struct();
    end

    settings = parse_settings( ...
        outerOptions, predictorOptions, FFprimeModel0.nBasis);

    context = struct();
    context.profilePsiN = profilePsiN;
    context.qLevels = problem.qEvaluation.interiorPsiN(:);
    context.contourNR = problem.qEvaluation.contourNR;
    context.fitRange = settings.fitRange;
    context.axisWeight = settings.axisWeight;
    context.minimumQPoints = settings.minimumQPoints;
    context.requiredAxisMode = '';
    context.requiredDpsiSign = 0;

    outerTimer = tic;
    evaluationCount = 0;

    % Select a valid cylindrical/user-supplied initial prediction. Scaling
    % candidates include FF'=0 as a safe fallback.
    initialCoefficients = FFprimeModel0.coefficients(:);
    initialScales = settings.initialScales(:).';

    if norm(initialCoefficients,inf) == 0
        initialScales = 1;
    else
        initialScales = unique(initialScales,'stable');
    end

    initialCandidates = cell(numel(initialScales),1);
    bestInitialIndex = 0;
    bestInitialMerit = inf;

    for k = 1:numel(initialScales)
        trialCoefficients = initialScales(k)*initialCoefficients;
        initialCandidates{k} = evaluate_candidate( ...
            problem,pprime,qTarget,FFprimeModel0,trialCoefficients, ...
            psi0,innerOptions,context);
        evaluationCount = evaluationCount+1;

        if initialCandidates{k}.valid && ...
                initialCandidates{k}.merit < bestInitialMerit
            bestInitialIndex = k;
            bestInitialMerit = initialCandidates{k}.merit;
        end
    end

    if bestInitialIndex == 0
        failureMessages = cellfun( ...
            @(candidate) candidate.failureReason, ...
            initialCandidates,'UniformOutput',false);
        error('solve_pq_constraint:NoValidInitialEquilibrium', ...
            ['No initial FFprime prediction produced a valid ', ...
             'equilibrium. Last failure: %s'],failureMessages{end});
    end

    current = initialCandidates{bestInitialIndex};
    currentModel = FFprimeModel0;
    currentModel.coefficients = current.coefficients;
    context.requiredAxisMode = current.state.axisMode;
    context.requiredDpsiSign = sign(current.state.dpsi);

    history = initialize_history(settings.maxIterations);
    acceptedIterations = 1;
    history = record_history(history,acceptedIterations,current, ...
        0,NaN,initialScales(bestInitialIndex),evaluationCount);

    converged = is_q_converged(current.metrics,settings) && ...
        acceptedIterations >= settings.minimumIterations;
    terminationReason = 'maximum_iterations';
    lastUpdate = empty_update();

    if settings.verbose
        print_iteration(acceptedIterations,current,0,NaN,evaluationCount);
    end

    % Actual-response reduced Gauss-Newton/LM iteration.
    while ~converged && acceptedIterations < settings.maxIterations
        [responseMatrix,differenceInfo,differenceEvaluations] = ...
            finite_difference_response( ...
                problem,pprime,qTarget,currentModel,current, ...
                innerOptions,context,settings);
        evaluationCount = evaluationCount+differenceEvaluations;

        if ~differenceInfo.success
            terminationReason = 'FFprime_response_failure';
            break
        end

        accepted = false;
        damping = settings.initialDamping;
        bestTrial = struct();
        bestTrialMerit = inf;
        acceptedStep = [];
        acceptedScale = NaN;
        acceptedPredictedReduction = NaN;

        for dampingAttempt = 1:settings.maximumDampingAttempts
            [coefficientStep,predictedReduction] = ...
                calculate_lm_step( ...
                    current,responseMatrix,damping,settings);

            if norm(coefficientStep,inf) > settings.maxCoefficientStep
                coefficientStep = coefficientStep ...
                    *(settings.maxCoefficientStep ...
                      /norm(coefficientStep,inf));
            end

            if norm(coefficientStep,inf) <= ...
                    settings.minimumCoefficientStep
                damping = damping*settings.dampingGrowth;
                continue
            end

            for lineIndex = 0:settings.maximumLineSearchSteps-1
                lineScale = 0.5^lineIndex;
                trialCoefficients = current.coefficients ...
                    +lineScale*coefficientStep;

                trial = evaluate_candidate( ...
                    problem,pprime,qTarget,currentModel, ...
                    trialCoefficients,current.state.psi, ...
                    innerOptions,context);
                evaluationCount = evaluationCount+1;

                if trial.valid && trial.merit < bestTrialMerit
                    bestTrial = trial;
                    bestTrialMerit = trial.merit;
                    acceptedStep = lineScale*coefficientStep;
                    acceptedScale = lineScale;
                    acceptedPredictedReduction = ...
                        lineScale*predictedReduction;
                end

                if trial.valid && is_acceptable_trial( ...
                        current,trial,lineScale, ...
                        predictedReduction,settings)
                    accepted = true;
                    break
                end
            end

            if accepted
                break
            end

            damping = damping*settings.dampingGrowth;
        end

        % Retain a strictly improving valid trial even if the approximate
        % response model missed the requested Armijo decrease.
        if ~accepted && ~isempty(fieldnames(bestTrial)) && ...
                bestTrialMerit < current.merit ...
                *(1-settings.minimumRelativeMeritDecrease)
            accepted = true;
        end

        if ~accepted
            terminationReason = 'no_acceptable_FFprime_step';
            break
        end

        previousMerit = current.merit;
        current = bestTrial;
        currentModel.coefficients = current.coefficients;
        acceptedIterations = acceptedIterations+1;

        lastUpdate = struct();
        lastUpdate.coefficientStep = acceptedStep;
        lastUpdate.lineSearchScale = acceptedScale;
        lastUpdate.predictedMeritReduction = ...
            acceptedPredictedReduction;
        lastUpdate.actualMeritReduction = previousMerit-current.merit;
        lastUpdate.responseMatrix = responseMatrix;
        lastUpdate.finiteDifference = differenceInfo;
        lastUpdate.damping = damping;

        history = record_history( ...
            history,acceptedIterations,current, ...
            norm(acceptedStep,inf),damping,acceptedScale,evaluationCount);

        if settings.verbose
            print_iteration(acceptedIterations,current, ...
                norm(acceptedStep,inf),acceptedScale,evaluationCount);
        end

        converged = is_q_converged(current.metrics,settings) && ...
            acceptedIterations >= settings.minimumIterations;

        if converged
            terminationReason = 'converged';
        end
    end

    history = trim_history(history,acceptedIterations);

    % Return one internally consistent last-accepted equilibrium.
    solution = current.state;
    solution.FFprimeModel = currentModel;
    solution.profileTable = current.profile;
    solution.qGeometry = current.qGeometry;
    solution.qResult = current.qResult;
    solution.qTarget = current.qTarget;
    solution.qRelativeResidual = current.metrics.relativeResidual;
    solution.qLogResidual = current.metrics.logResidual;
    solution.qValidMask = current.metrics.validMask;

    result = struct();
    result.converged = converged;
    result.terminationReason = terminationReason;
    result.outerIterations = acceptedIterations;
    result.functionEvaluations = evaluationCount;
    result.elapsedTime = toc(outerTimer);
    result.finalMaxAbsRelativeQ = current.metrics.maxAbsRelative;
    result.finalRmsLogQ = current.metrics.rmsLog;
    result.finalValidQPointCount = current.metrics.validPointCount;
    result.finalFFprimeModel = currentModel;
    result.finalFProfile = current.profile;
    result.finalInnerResult = current.innerResult;
    result.lastAppliedFFprimeUpdate = lastUpdate;
    result.history = history;
    result.tolerances = struct( ...
        'qRelative',settings.qRelativeTolerance, ...
        'qRmsLog',settings.qRmsLogTolerance);

    if ~converged && settings.errorOnFailure
        error('solve_pq_constraint:Nonconvergence', ...
            ['The p-q solve did not converge after %d accepted ', ...
             'iterations. Termination reason: %s.'], ...
            acceptedIterations,terminationReason);
    end
end


function candidate = evaluate_candidate( ...
    problem,pprime,qTarget,model,coefficients,psiGuess, ...
    innerOptions,context)
%EVALUATE_CANDIDATE Perform one complete FF' -> GS -> q evaluation.

    candidate = empty_candidate(coefficients);
    trialModel = model;
    trialModel.coefficients = coefficients(:);

    try
        profile = evaluate_FFprime_model( ...
            trialModel,context.profilePsiN);
        profiles = struct( ...
            'pprime',pprime(:),'FFprime',profile.FFprime(:));

        warningState = warning('query','GS:PicardPrepared:NoConvergence');
        warning('off','GS:PicardPrepared:NoConvergence');
        warningCleanup = onCleanup(@() warning( ...
            warningState.state,'GS:PicardPrepared:NoConvergence'));
        [state,innerResult] = solve_GS_inner( ...
            problem,profiles,psiGuess,innerOptions);
        clear warningCleanup

        if ~innerResult.converged
            candidate.failureReason = 'inner_GS_nonconvergence';
            candidate.state = state;
            candidate.innerResult = innerResult;
            return
        end

        if ~isempty(context.requiredAxisMode) && ...
                ~strcmp(state.axisMode,context.requiredAxisMode)
            candidate.failureReason = 'magnetic_axis_branch_change';
            return
        end

        if context.requiredDpsiSign ~= 0 && ...
                sign(state.dpsi) ~= context.requiredDpsiSign
            candidate.failureReason = 'flux_span_sign_change';
            return
        end

        qGeometry = evaluate_q_geometry( ...
            problem.mesh,state.psi,state.psiAxis,state.dpsi, ...
            context.qLevels,context.contourNR, ...
            state.axisNode,state.axisMode);

        warningState = warning('query', ...
            'GS:QBoundary:ExtrapolationSensitive');
        warning('off','GS:QBoundary:ExtrapolationSensitive');
        warningCleanup = onCleanup(@() warning( ...
            warningState.state,'GS:QBoundary:ExtrapolationSensitive'));
        qResult = evaluate_q_from_G( ...
            qGeometry,profile.psiN,profile.G,profile.signF);
        clear warningCleanup

        qPsiN = qGeometry.psiN(:);
        qValues = qResult.q(:);
        qGeometryFactor = qGeometry.qPerF(:);
        qTargetValues = evaluate_q_target(qTarget,qPsiN);

        metrics = calculate_q_metrics( ...
            qPsiN,qValues,qTargetValues,qGeometryFactor, ...
            context.fitRange,context.axisWeight, ...
            context.minimumQPoints);

        candidate.valid = true;
        candidate.failureReason = '';
        candidate.coefficients = coefficients(:);
        candidate.profile = profile;
        candidate.state = state;
        candidate.innerResult = innerResult;
        candidate.qGeometry = qGeometry;
        candidate.qResult = qResult;
        candidate.qTarget = qTargetValues;
        candidate.metrics = metrics;
        candidate.residual = metrics.logResidual(metrics.validMask);
        candidate.weights = metrics.weights(metrics.validMask);
        candidate.merit = metrics.weightedMerit;
    catch ME
        candidate.failureReason = sprintf('%s: %s', ...
            ME.identifier,ME.message);
    end
end


function [response,info,nEvaluation] = finite_difference_response( ...
    problem,pprime,qTarget,model,current,innerOptions,context,settings)
%FINITE_DIFFERENCE_RESPONSE Actual reduced q response to FF' coefficients.

    nCoefficient = numel(current.coefficients);
    nResidual = numel(current.residual);
    response = nan(nResidual,nCoefficient);
    stepUsed = nan(nCoefficient,1);
    directionUsed = zeros(nCoefficient,1);
    nEvaluation = 0;

    for j = 1:nCoefficient
        step = max( ...
            settings.finiteDifferenceAbsoluteStep, ...
            settings.finiteDifferenceRelativeStep ...
                *max(1,abs(current.coefficients(j))));
        success = false;

        for direction = [1,-1]
            trialCoefficients = current.coefficients;
            trialCoefficients(j) = ...
                trialCoefficients(j)+direction*step;
            trial = evaluate_candidate( ...
                problem,pprime,qTarget,model,trialCoefficients, ...
                current.state.psi,innerOptions,context);
            nEvaluation = nEvaluation+1;

            if trial.valid && numel(trial.residual) == nResidual && ...
                    isequal(trial.metrics.validMask, ...
                            current.metrics.validMask)
                response(:,j) = ...
                    (trial.residual-current.residual)/(direction*step);
                stepUsed(j) = step;
                directionUsed(j) = direction;
                success = true;
                break
            end
        end

        if ~success
            info = struct( ...
                'success',false,'failedCoefficient',j, ...
                'stepUsed',stepUsed,'directionUsed',directionUsed);
            return
        end
    end

    info = struct();
    info.success = all(isfinite(response(:)));
    info.failedCoefficient = NaN;
    info.stepUsed = stepUsed;
    info.directionUsed = directionUsed;
    info.rcond = rcond(response.'*response ...
        +settings.jacobianRidge*eye(nCoefficient));
end


function [step,predictedReduction] = calculate_lm_step( ...
    current,response,damping,settings)
%CALCULATE_LM_STEP Weighted regularized Levenberg-Marquardt step.

    sqrtWeight = sqrt(current.weights(:));
    weightedResponse = response.*sqrtWeight;
    weightedResidual = current.residual.*sqrtWeight;
    nCoefficient = numel(current.coefficients);
    regularizationMatrix = diag((0:nCoefficient-1).');
    gaussNewtonMatrix = weightedResponse.'*weightedResponse;
    curvatureScale = max( ...
        trace(gaussNewtonMatrix)/max(nCoefficient,1),1.0e-10);
    dampingMatrix = diag(max( ...
        diag(gaussNewtonMatrix),curvatureScale));

    normalMatrix = gaussNewtonMatrix ...
        +settings.regularization ...
            *(regularizationMatrix.'*regularizationMatrix) ...
        +damping*dampingMatrix ...
        +settings.jacobianRidge*eye(nCoefficient);

    gradient = weightedResponse.'*weightedResidual ...
        +settings.regularization ...
            *(regularizationMatrix.'*regularizationMatrix) ...
            *current.coefficients;

    step = -normalMatrix\gradient;
    predictedResidual = weightedResidual+weightedResponse*step;
    predictedReduction = 0.5*( ...
        weightedResidual.'*weightedResidual ...
        -predictedResidual.'*predictedResidual) ...
        /max(sum(current.weights),eps);
end


function accepted = is_acceptable_trial( ...
    current,trial,lineScale,predictedReduction,settings)
%IS_ACCEPTABLE_TRIAL Actual-GS sufficient-decrease test.

    requestedDecrease = max( ...
        settings.minimumAbsoluteMeritDecrease, ...
        settings.armijoFraction*lineScale*max(predictedReduction,0));
    accepted = trial.merit <= current.merit-requestedDecrease;

    if trial.metrics.maxAbsRelative <= settings.qRelativeTolerance && ...
            trial.metrics.rmsLog <= settings.qRmsLogTolerance
        accepted = true;
    end
end


function metrics = calculate_q_metrics( ...
    psiN,qValues,qTarget,qPerF,fitRange,axisWeight,minimumQPoints)
%CALCULATE_Q_METRICS Construct optimization and convergence metrics.

    psiN = psiN(:);
    qValues = qValues(:);
    qTarget = qTarget(:);
    qPerF = qPerF(:);

    validMask = isfinite(psiN) & isfinite(qValues) ...
        &isfinite(qTarget) & isfinite(qPerF) ...
        &qValues ~= 0 & qTarget ~= 0 & qPerF ~= 0 ...
        &psiN >= fitRange(1) & psiN <= fitRange(2);

    if any(validMask & sign(qValues) ~= sign(qTarget))
        error('solve_pq_constraint:InconsistentQSign', ...
            'Calculated and prescribed q have inconsistent signs.');
    end

    if nnz(validMask) < minimumQPoints
        error('solve_pq_constraint:TooFewValidQPoints', ...
            'Too few valid q surfaces are available.');
    end

    relativeResidual = nan(size(qValues));
    logResidual = nan(size(qValues));
    weights = zeros(size(qValues));
    relativeResidual(validMask) = ...
        (qValues(validMask)-qTarget(validMask))./qTarget(validMask);
    logResidual(validMask) = ...
        log(abs(qValues(validMask)./qTarget(validMask)));
    weights(validMask) = 1;
    axisMask = validMask & abs(psiN) <= 100*eps;
    weights(axisMask) = axisWeight;

    validWeights = weights(validMask);
    validLogResidual = logResidual(validMask);
    metrics = struct();
    metrics.validMask = validMask;
    metrics.validPointCount = nnz(validMask);
    metrics.relativeResidual = relativeResidual;
    metrics.logResidual = logResidual;
    metrics.weights = weights;
    metrics.maxAbsRelative = ...
        max(abs(relativeResidual(validMask)));
    metrics.rmsLog = sqrt(mean(validLogResidual.^2));
    metrics.weightedMerit = 0.5*sum( ...
        validWeights.*validLogResidual.^2)/sum(validWeights);
end


function values = evaluate_q_target(qTarget,psiN)
%EVALUATE_Q_TARGET Evaluate a supported prescribed-q representation.

    if isa(qTarget,'function_handle')
        values = qTarget(psiN);
    elseif isstruct(qTarget)
        if ~isfield(qTarget,'psiN') || ~isfield(qTarget,'q')
            error('solve_pq_constraint:QTargetStructure', ...
                'qTarget must contain psiN and q.');
        end
        values = interp1( ...
            qTarget.psiN(:),qTarget.q(:),psiN,'pchip');
    elseif isscalar(qTarget)
        values = qTarget*ones(size(psiN));
    elseif isnumeric(qTarget) && numel(qTarget) == numel(psiN)
        values = qTarget(:);
    else
        error('solve_pq_constraint:QTarget', ...
            'Unsupported qTarget representation.');
    end

    values = values(:);

    if numel(values) ~= numel(psiN) || ...
            any(~isfinite(values)) || any(values == 0)
        error('solve_pq_constraint:InvalidQTarget', ...
            'The evaluated q target must be finite and nonzero.');
    end
end


function settings = parse_settings(outer,predictor,nBasis)
%PARSE_SETTINGS Validate reduced inverse-solver controls.

    settings = struct();
    settings.maxIterations = get_option(outer,'maxIterations',20);
    settings.minimumIterations = get_option(outer,'minimumIterations',1);
    settings.qRelativeTolerance = ...
        get_option(outer,'qRelativeTolerance',1.0e-3);
    settings.qRmsLogTolerance = ...
        get_option(outer,'qRmsLogTolerance',1.0e-3);
    settings.minimumQPoints = get_option(outer,'minimumQPoints',3);
    settings.errorOnFailure = logical( ...
        get_option(outer,'errorOnFailure',false));
    settings.verbose = logical(get_option(outer,'verbose',true));

    settings.fitRange = get_option(predictor,'fitRange',[0,0.95]);
    settings.axisWeight = get_option(predictor,'axisWeight',0.25);
    settings.initialScales = ...
        get_option(predictor,'initialScales',[1,0.5,0]);
    settings.finiteDifferenceAbsoluteStep = ...
        get_option(predictor,'finiteDifferenceAbsoluteStep',0.05);
    settings.finiteDifferenceRelativeStep = ...
        get_option(predictor,'finiteDifferenceRelativeStep',0.02);
    settings.regularization = ...
        get_option(predictor,'regularization',1.0e-4);
    settings.jacobianRidge = ...
        get_option(predictor,'jacobianRidge',1.0e-10);
    settings.initialDamping = ...
        get_option(predictor,'initialDamping',1.0e-2);
    settings.dampingGrowth = get_option(predictor,'dampingGrowth',10);
    settings.maximumDampingAttempts = ...
        get_option(predictor,'maximumDampingAttempts',4);
    settings.maximumLineSearchSteps = ...
        get_option(predictor,'maximumLineSearchSteps',6);
    settings.maxCoefficientStep = ...
        get_option(predictor,'maxCoefficientStep',2.0);
    settings.minimumCoefficientStep = ...
        get_option(predictor,'minimumCoefficientStep',1.0e-8);
    settings.armijoFraction = ...
        get_option(predictor,'armijoFraction',1.0e-3);
    settings.minimumAbsoluteMeritDecrease = ...
        get_option(predictor,'minimumAbsoluteMeritDecrease',1.0e-12);
    settings.minimumRelativeMeritDecrease = ...
        get_option(predictor,'minimumRelativeMeritDecrease',1.0e-8);

    validateattributes(settings.maxIterations,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(settings.minimumIterations,{'numeric'}, ...
        {'scalar','integer','positive','<=',settings.maxIterations});
    validateattributes(settings.qRelativeTolerance,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(settings.qRmsLogTolerance,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(settings.minimumQPoints,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(settings.axisWeight,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(settings.finiteDifferenceAbsoluteStep,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(settings.finiteDifferenceRelativeStep,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(settings.regularization,{'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    validateattributes(settings.jacobianRidge,{'numeric'}, ...
        {'scalar','real','finite','nonnegative'});
    validateattributes(settings.initialDamping,{'numeric'}, ...
        {'scalar','real','finite','positive'});
    validateattributes(settings.dampingGrowth,{'numeric'}, ...
        {'scalar','real','finite','>',1});
    validateattributes(settings.maximumDampingAttempts,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(settings.maximumLineSearchSteps,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(settings.maxCoefficientStep,{'numeric'}, ...
        {'scalar','real','finite','positive'});

    if numel(settings.fitRange) ~= 2 || ...
            any(~isfinite(settings.fitRange)) || ...
            settings.fitRange(1) < 0 || settings.fitRange(2) > 1 || ...
            settings.fitRange(1) >= settings.fitRange(2)
        error('solve_pq_constraint:FitRange', ...
            'FFprime.fitRange must be an increasing interval in [0,1].');
    end

    if isempty(settings.initialScales) || ...
            any(~isfinite(settings.initialScales)) || ...
            any(settings.initialScales < 0)
        error('solve_pq_constraint:InitialScales', ...
            'FFprime.initialScales must be finite and nonnegative.');
    end

    if nBasis < 1
        error('solve_pq_constraint:BasisCount', ...
            'At least one FFprime basis coefficient is required.');
    end
end


function history = initialize_history(maxIterations)
%INITIALIZE_HISTORY Allocate accepted-iterate diagnostics.

    history = struct();
    history.iteration = nan(maxIterations,1);
    history.innerIterations = nan(maxIterations,1);
    history.innerResidual = nan(maxIterations,1);
    history.qMaxRelativeError = nan(maxIterations,1);
    history.qRmsLogError = nan(maxIterations,1);
    history.merit = nan(maxIterations,1);
    history.coefficientStep = nan(maxIterations,1);
    history.damping = nan(maxIterations,1);
    history.lineSearchScale = nan(maxIterations,1);
    history.functionEvaluations = nan(maxIterations,1);
    history.coefficients = cell(maxIterations,1);
    history.FFprime = cell(maxIterations,1);
end


function history = record_history( ...
    history,index,candidate,coefficientStep,damping,lineScale,nEvaluation)
%RECORD_HISTORY Store one accepted reduced iterate.

    history.iteration(index) = index;
    history.innerIterations(index) = candidate.innerResult.iterations;
    history.innerResidual(index) = ...
        candidate.innerResult.finalRelativeResidual;
    history.qMaxRelativeError(index) = ...
        candidate.metrics.maxAbsRelative;
    history.qRmsLogError(index) = candidate.metrics.rmsLog;
    history.merit(index) = candidate.merit;
    history.coefficientStep(index) = coefficientStep;
    history.damping(index) = damping;
    history.lineSearchScale(index) = lineScale;
    history.functionEvaluations(index) = nEvaluation;
    history.coefficients{index} = candidate.coefficients;
    history.FFprime{index} = candidate.profile.FFprime;
end


function history = trim_history(history,nIteration)
%TRIM_HISTORY Remove unused preallocation rows.

    fieldNames = fieldnames(history);

    for k = 1:numel(fieldNames)
        values = history.(fieldNames{k});
        history.(fieldNames{k}) = values(1:nIteration,:);
    end
end


function converged = is_q_converged(metrics,settings)
%IS_Q_CONVERGED Apply both prescribed-q tolerances.

    converged = ...
        metrics.maxAbsRelative <= settings.qRelativeTolerance && ...
        metrics.rmsLog <= settings.qRmsLogTolerance;
end


function print_iteration(iteration,candidate,step,lineScale,nEvaluation)
%PRINT_ITERATION Concise actual-GS outer diagnostics.

    fprintf([ ...
        '[p-q FFprime] %3d: inner %3d, max|dq/q| %.3e, ', ...
        'rms(log q/q*) %.3e, merit %.3e, ', ...
        '|dc|inf %.3e, alpha %.3g, eval %d\n'], ...
        iteration,candidate.innerResult.iterations, ...
        candidate.metrics.maxAbsRelative,candidate.metrics.rmsLog, ...
        candidate.merit,step,lineScale,nEvaluation);
end


function candidate = empty_candidate(coefficients)
%EMPTY_CANDIDATE Default invalid trial result.

    candidate = struct();
    candidate.valid = false;
    candidate.failureReason = 'not_evaluated';
    candidate.coefficients = coefficients(:);
    candidate.profile = struct();
    candidate.state = struct();
    candidate.innerResult = struct();
    candidate.qGeometry = struct();
    candidate.qResult = struct();
    candidate.qTarget = [];
    candidate.metrics = struct();
    candidate.residual = [];
    candidate.weights = [];
    candidate.merit = inf;
end


function update = empty_update()
%EMPTY_UPDATE Default before an accepted coefficient update exists.

    update = struct();
    update.coefficientStep = [];
    update.lineSearchScale = NaN;
    update.predictedMeritReduction = NaN;
    update.actualMeritReduction = NaN;
    update.responseMatrix = [];
    update.finiteDifference = struct();
    update.damping = NaN;
end


function substructure = get_substructure(value,fieldName)
%GET_SUBSTRUCTURE Read and validate an optional substructure.

    if isfield(value,fieldName) && ~isempty(value.(fieldName))
        substructure = value.(fieldName);
    else
        substructure = struct();
    end

    if ~isstruct(substructure) || ~isscalar(substructure)
        error('solve_pq_constraint:InvalidOptions', ...
            'options.%s must be a scalar structure.',fieldName);
    end
end


function value = get_option(options,fieldName,defaultValue)
%GET_OPTION Read an optional field.

    if isfield(options,fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
    else
        value = defaultValue;
    end
end


function assert_required_fields(value,requiredFields,valueName)
%ASSERT_REQUIRED_FIELDS Validate a scalar structure interface.

    if ~isstruct(value) || ~isscalar(value)
        error('solve_pq_constraint:InvalidStructure', ...
            '%s must be a scalar structure.',valueName);
    end

    for k = 1:numel(requiredFields)
        if ~isfield(value,requiredFields{k})
            error('solve_pq_constraint:MissingField', ...
                '%s.%s is required.',valueName,requiredFields{k});
        end
    end
end
