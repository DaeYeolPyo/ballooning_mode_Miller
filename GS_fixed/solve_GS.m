function equilibrium = solve_GS(input)
    %======================================================
    % Set dependencies
    %======================================================
    setup_solve_GS;
    
    %======================================================
    % Prepare profile-independent GS problem
    %======================================================
    problem = prepare_GS_problem(input);

    %======================================================
    % Load prepared objects
    %======================================================
    P3 = problem.mesh;

    R0 = problem.dim.R0;
    B0 = problem.dim.B0;

    psiNProfile = problem.profile.psiN;

    verbose = problem.output.verbose;
    checkTime = problem.output.checkTime;
    showPlot = problem.output.showPlot;

    %======================================================
    % Profile inputs and solver dispatch
    %======================================================
    if ~isfield(input, 'profile') || ...
            ~isfield(input.profile, 'mode')
        error('GS:Solve:MissingProfileMode', ...
            'input.profile.mode is required.');
    end

    profileMode = input.profile.mode;

    validateattributes(profileMode, {'numeric'}, ...
        {'scalar', 'integer'});

    if ~ismember(profileMode, [1, 3])
        error('GS:Solve:UnsupportedProfileMode', ...
            ['Supported profile modes are currently ', ...
             '1 (pprime-FFprime) and 3 (pprime-q).']);
    end

    if ~isfield(input.profile, 'pprime')
        error('GS:Solve:MissingPprime', ...
            'input.profile.pprime is required.');
    end

    if ~isfield(input.profile, 'pedge')
        error('GS:Solve:MissingPedge', ...
            'input.profile.pedge is required.');
    end

    pprime = input.profile.pprime(:);
    pedge = input.profile.pedge;

    if numel(pprime) ~= numel(psiNProfile)
        error('GS:Solve:PprimeSize', ...
            ['input.profile.pprime must have the same length as ', ...
             'problem.profile.psiN.']);
    end

    constraintResult = [];
    FFprimeModel = [];

    solverTimer = tic;

    switch profileMode
        %==============================================================
        % Mode 1: prescribed pprime and FFprime
        %==============================================================
        case 1
            if ~isfield(input.profile, 'FFprime')
                error('GS:Solve:MissingFFprime', ...
                    ['input.profile.FFprime is required for ', ...
                     'profile mode 1.']);
            end

            FFprime = input.profile.FFprime(:);

            if numel(FFprime) ~= numel(psiNProfile)
                error('GS:Solve:FFprimeSize', ...
                    ['input.profile.FFprime must have the same ', ...
                     'length as problem.profile.psiN.']);
            end

            profiles = struct();
            profiles.pprime = pprime;
            profiles.FFprime = FFprime;

            [state, PicardResult] = solve_GS_inner( ...
                problem, ...
                profiles, ...
                problem.initial.psi, ...
                problem.options.Picard);

            % Recover G and F from the prescribed FFprime.
            Fvac = R0*B0;

            FProfile = reconstruct_G_from_FFprime( ...
                psiNProfile, ...
                FFprime, ...
                Fvac);

            % Evaluate q for the converged mode-1 equilibrium.
            qGeometry = evaluate_q_geometry( ...
                P3, ...
                state.psi, ...
                state.psiAxis, ...
                state.dpsi, ...
                problem.qEvaluation.interiorPsiN, ...
                problem.qEvaluation.contourNR, ...
                state.axisNode, ...
                state.axisMode);

            qResult = evaluate_q_from_G( ...
                qGeometry, ...
                FProfile.psiN, ...
                FProfile.G, ...
                FProfile.signF);

        %==============================================================
        % Mode 3: prescribed pprime and q
        %==============================================================
        case 3
            if ~isfield(input.profile, 'q')
                error('GS:Solve:MissingQ', ...
                    ['input.profile.q is required for ', ...
                     'profile mode 3.']);
            end

            % Convert the solve_GS input format into a qTarget format
            % accepted by solve_pq_constraint.
            qTarget = normalize_q_target_input( ...
                input.profile.q, psiNProfile);

            pqOptions = make_pq_constraint_options( ...
                input, problem, verbose);

            % Construct the cylindrical/direct-FFprime initial model.
            % initialize_FFprime_model also reads optional fields:
            %   input.profile.normalization
            %   input.profile.FSign
            %   input.profile.initial
            FFprimeModel0 = initialize_FFprime_model( ...
                problem, input.profile, pprime, qTarget, ...
                pqOptions.FFprime);

            [state, constraintResult] = solve_pq_constraint( ...
                problem, ...
                pprime, ...
                qTarget, ...
                FFprimeModel0, ...
                problem.initial.psi, ...
                pqOptions);

            % Keep the final inner Picard result as solverResult so that
            % existing report and plotting functions remain compatible.
            PicardResult = ...
                constraintResult.finalInnerResult;

            FFprimeModel = ...
                constraintResult.finalFFprimeModel;

            FProfile = constraintResult.finalFProfile;

            % q has already been evaluated for the final outer iterate.
            qGeometry = state.qGeometry;
            qResult = state.qResult;

            if ~isstruct(qResult) || ...
                    ~isfield(qResult, 'q') || ...
                    isempty(qResult.q)
                error('GS:Solve:MissingFinalQ', ...
                    ['The p-q solver terminated before producing ', ...
                     'a valid q profile.']);
            end
    end

    solveTime = toc(solverTimer);

    if checkTime
        fprintf('Solver time           = %.3f s\n', solveTime);
    end

    %======================================================
    % Solver convergence checks
    %======================================================
    if ~PicardResult.converged
        warning('GS:Solve:InnerNonconvergence', ...
            ['The final inner Grad-Shafranov solve did not ', ...
             'converge.']);
    end

    if profileMode == 3 && ~constraintResult.converged
        warning('GS:Solve:PQNonconvergence', ...
            ['The p-q constraint iteration did not converge. ', ...
             'Termination reason: %s.'], ...
            constraintResult.terminationReason);
    end

    %======================================================
    % Normalize mode-dependent results
    %======================================================
    solved = struct();

    solved.mode = profileMode;

    solved.state = state;
    solved.innerResult = PicardResult;

    solved.pprime = pprime;
    solved.pedge = pedge;

    solved.FProfile = FProfile;

    solved.qGeometry = qGeometry;
    solved.qResult = qResult;

    solved.solveTime = solveTime;

    %======================================================
    % Mode-3-specific results
    %======================================================
    if profileMode == 3
        solved.constraintResult = constraintResult;
        solved.FFprimeModel = FFprimeModel;

        if ~isfield(state, 'qTarget') || isempty(state.qTarget)
            error('GS:Solve:MissingFinalQTarget', ...
                ['The p-q solver did not return the prescribed q ', ...
                 'values on its final q-evaluation grid.']);
        end

        solved.qTargetEvaluated = state.qTarget;
    end

    %======================================================
    % Common equilibrium finalization
    %======================================================
    equilibrium = finalize_GS_equilibrium( ...
        problem, solved);

    %======================================================
    % Report and plot
    %=====================================================
    if verbose
        report_GS_result(equilibrium);
    end

    if showPlot
        plot_GS_equilibrium(equilibrium);
    end
end

function qTarget = normalize_q_target_input( ...
    qInput, psiNProfile)
%NORMALIZE_Q_TARGET_INPUT
% Convert solve_GS profile input into the representation expected by
% solve_pq_constraint.

    if isa(qInput, 'function_handle')
        qTarget = qInput;
        return
    end

    if isstruct(qInput)
        qTarget = qInput;
        return
    end

    if ~isnumeric(qInput)
        error('GS:Solve:InvalidQInput', ...
            ['input.profile.q must be numeric, a function handle, ', ...
             'or a structure containing psiN and q.']);
    end

    if isscalar(qInput)
        qTarget = qInput;
        return
    end

    qInput = qInput(:);

    if numel(qInput) ~= numel(psiNProfile)
        error('GS:Solve:QProfileSize', ...
            ['A numeric input.profile.q vector must have the ', ...
             'same length as problem.profile.psiN.']);
    end

    qTarget = struct();
    qTarget.psiN = psiNProfile(:);
    qTarget.q = qInput;
end


function options = make_pq_constraint_options( ...
    input, problem, verbose)
%MAKE_PQ_CONSTRAINT_OPTIONS
% Construct options for solve_pq_constraint.

    if isfield(input, 'pqConstraint') && ...
            ~isempty(input.pqConstraint)

        if ~isstruct(input.pqConstraint) || ...
                ~isscalar(input.pqConstraint)
            error('GS:Solve:InvalidPQOptions', ...
                'input.pqConstraint must be a scalar structure.');
        end

        options = input.pqConstraint;
    else
        options = struct();
    end

    if ~isfield(options, 'outer') || ...
            isempty(options.outer)
        options.outer = struct();
    end

    if ~isstruct(options.outer) || ...
            ~isscalar(options.outer)
        error('GS:Solve:InvalidPQOuterOptions', ...
            'input.pqConstraint.outer must be a scalar structure.');
    end

    if ~isfield(options.outer, 'verbose') || ...
            isempty(options.outer.verbose)
        options.outer.verbose = verbose;
    end

    % input.Picard remains the canonical source of inner-solver options,
    % but trial solves follow the outer-loop verbosity setting.
    options.inner = problem.options.Picard;
    options.inner.verbose = logical(options.outer.verbose);

    if isfield(options, 'updateG')
        error('GS:Solve:ObsoleteGUpdateOptions', ...
            ['input.pqConstraint.updateG belonged to the removed ', ...
             'frozen-logG algorithm. Use input.pqConstraint.FFprime.']);
    end

    if ~isfield(options, 'FFprime') || isempty(options.FFprime)
        options.FFprime = struct();
    end

    if ~isstruct(options.FFprime) || ~isscalar(options.FFprime)
        error('GS:Solve:InvalidFFprimeOptions', ...
            'input.pqConstraint.FFprime must be a scalar structure.');
    end
end
