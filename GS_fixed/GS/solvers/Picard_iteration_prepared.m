function [psi, result] = Picard_iteration_prepared( ...
    problem, pprime, FFprime, psi0, options)
%PICARD_ITERATION_PREPARED
%
% Solve the nonlinear fixed-boundary Grad-Shafranov equation using
% preassembled and prefactorized FEM operators.
%
% The following objects are reused from prepare_GS_problem:
%
%   problem.mesh
%   problem.quad
%   problem.operator.K
%   problem.operator.Kfactor
%   problem.operator.freeNodes
%   problem.operator.boundaryNodes
%   problem.operator.boundaryValues
%   problem.operator.boundaryCorrection
%
%
% INPUT
%
%   problem
%       Output of prepare_GS_problem.
%
%   pprime
%       dp/dpsi_N samples on the uniform normalized-flux grid.
%
%   FFprime
%       F*dF/dpsi_N samples on the same grid.
%
%   psi0
%       Initial dimensional poloidal-flux vector.
%
%   options
%       Picard iteration options. If omitted, use
%
%           problem.options.Picard.
%
%
% OUTPUT
%
%   psi
%       Converged dimensional poloidal flux.
%
%   result
%       Picard convergence information.

    narginchk(4, 5);

    validate_prepared_problem(problem);

    if nargin < 5 || isempty(options)
        if isfield(problem, 'options') && ...
                isfield(problem.options, 'Picard')
            options = problem.options.Picard;
        else
            options = struct();
        end
    end

    opts = parse_Picard_options(options);

    %==============================================================
    % Prepared objects
    %==============================================================
    P3 = problem.mesh;
    quad = problem.quad;

    K = problem.operator.K;
    Kfactor = problem.operator.Kfactor;

    freeNodes = ...
        problem.operator.freeNodes(:);

    boundaryNodes = ...
        problem.operator.boundaryNodes(:);

    psiB = ...
        problem.operator.boundaryValues(:);

    boundaryCorrection = ...
        problem.operator.boundaryCorrection(:);

    symmetryError = ...
        problem.operator.symmetryError;

    psiBoundary = ...
        problem.boundary.psi;

    Ndof = size(P3.points,1);

    %==============================================================
    % Current RHS implementation supports psi_boundary = 0 only
    %==============================================================
    boundaryTolerance = ...
        100*eps(max(1,abs(psiBoundary)));

    if abs(psiBoundary) > boundaryTolerance
        error('GS:PicardPrepared:NonzeroBoundary', ...
            ['construct_RHS currently assumes psi_boundary = 0. ', ...
             'Generalize construct_RHS before using a nonzero ', ...
             'boundary flux.']);
    end

    %==============================================================
    % Profile validation
    %==============================================================
    pprime = pprime(:);
    FFprime = FFprime(:);

    if numel(pprime) ~= numel(FFprime)
        error('GS:PicardPrepared:ProfileSize', ...
            'pprime and FFprime must have equal lengths.');
    end

    if numel(pprime) < 2
        error('GS:PicardPrepared:ProfileSize', ...
            'At least two profile samples are required.');
    end

    if any(~isfinite(pprime))
        error('GS:PicardPrepared:InvalidPprime', ...
            'pprime contains NaN or Inf.');
    end

    if any(~isfinite(FFprime))
        error('GS:PicardPrepared:InvalidFFprime', ...
            'FFprime contains NaN or Inf.');
    end

    if isfield(problem, 'profile') && ...
            isfield(problem.profile, 'psiN')

        psiNProfile = problem.profile.psiN(:);

        if numel(psiNProfile) ~= numel(pprime)
            error('GS:PicardPrepared:ProfileGridSize', ...
                ['The profile length does not match ', ...
                 'problem.profile.psiN.']);
        end

        % construct_RHS currently creates its interpolation grid with
        % linspace(0,1,nProfile), so the prepared grid must be uniform.
        expectedPsiN = ...
            linspace(0,1,numel(psiNProfile)).';

        gridTolerance = ...
            100*eps(max(1,numel(psiNProfile)));

        if norm(psiNProfile-expectedPsiN,inf) > ...
                gridTolerance
            error('GS:PicardPrepared:NonuniformProfileGrid', ...
                ['construct_RHS currently requires a uniform ', ...
                 'normalized-flux profile grid.']);
        end
    end

    %==============================================================
    % Initial psi
    %==============================================================
    psi = psi0(:);

    if numel(psi) ~= Ndof
        error('GS:PicardPrepared:InitialPsiSize', ...
            'psi0 must contain one value per P3 degree of freedom.');
    end

    if any(~isfinite(psi))
        error('GS:PicardPrepared:InvalidInitialPsi', ...
            'psi0 contains NaN or Inf.');
    end

    if numel(psiB) ~= numel(boundaryNodes)
        error('GS:PicardPrepared:BoundarySize', ...
            ['The number of boundary values does not match ', ...
             'the number of boundary nodes.']);
    end

    initialBoundaryError = norm( ...
        psi(boundaryNodes)-psiB, inf);

    psi(boundaryNodes) = psiB;

    %==============================================================
    % Verify prepared operator compatibility
    %==============================================================
    if ~isequal(size(K), [Ndof,Ndof])
        error('GS:PicardPrepared:OperatorSize', ...
            'The prepared stiffness matrix has incompatible size.');
    end

    if numel(boundaryCorrection) ~= numel(freeNodes)
        error('GS:PicardPrepared:BoundaryCorrectionSize', ...
            'The prepared boundary correction has incompatible size.');
    end

    if symmetryError > opts.symmetryTolerance
        error('GS:PicardPrepared:Symmetry', ...
            ['Prepared KII symmetry error %.3e exceeds ', ...
             'the requested tolerance %.3e.'], ...
            symmetryError, opts.symmetryTolerance);
    end

    %==============================================================
    % Initial nonlinear source
    %==============================================================
    [psiAxis, dpsi, resolvedAxisMode] = ...
        evaluate_axis_and_span( ...
            psi, freeNodes, psiBoundary, ...
            opts.axisMode, ...
            opts.fluxSpanTolerance);

    fCurrent = assemble_RHS( ...
        P3, psi, quad, ...
        pprime, FFprime, dpsi);

    %==============================================================
    % History arrays
    %==============================================================
    maxIter = opts.maxIterations;

    updateHistory   = nan(maxIter,1);
    residualHistory = nan(maxIter,1);
    psiAxisHistory  = nan(maxIter,1);
    dpsiHistory     = nan(maxIter,1);
    minPsiNHistory  = nan(maxIter,1);
    maxPsiNHistory  = nan(maxIter,1);

    converged = false;

    residualVector = ...
        nan(numel(freeNodes),1);

    if opts.verbose
        fprintf('\n');
        fprintf('=======================================================\n');
        fprintf(' Prepared Grad-Shafranov Picard iteration\n');
        fprintf('=======================================================\n');
        fprintf([' iter     psi_axis     dpsi     ', ...
                 'update error     residual error\n']);
    end

    %==============================================================
    % Picard iteration
    %==============================================================
    for iter = 1:maxIter

        %----------------------------------------------------------
        % Solve with the lagged nonlinear source
        %----------------------------------------------------------
        rhsI = ...
            fCurrent(freeNodes)-boundaryCorrection;

        psiSolved = psi;

        psiSolved(freeNodes) = ...
            Kfactor\rhsI;

        psiSolved(boundaryNodes) = ...
            psiB;

        %----------------------------------------------------------
        % Under-relaxation
        %----------------------------------------------------------
        psiNew = psi;

        psiNew(freeNodes) = ...
            (1-opts.omega)*psi(freeNodes) ...
            +opts.omega*psiSolved(freeNodes);

        psiNew(boundaryNodes) = psiB;

        %----------------------------------------------------------
        % Relative solution update
        %----------------------------------------------------------
        updateVector = ...
            psiNew(freeNodes)-psi(freeNodes);

        updateScale = max([ ...
            norm(psiNew(freeNodes)), ...
            norm(psi(freeNodes)), ...
            eps]);

        relativeUpdate = ...
            norm(updateVector)/updateScale;

        %----------------------------------------------------------
        % Updated axis and signed flux span
        %----------------------------------------------------------
        [psiAxisNew, dpsiNew, resolvedAxisModeNew] = ...
            evaluate_axis_and_span( ...
                psiNew, freeNodes, psiBoundary, ...
                opts.axisMode, ...
                opts.fluxSpanTolerance);

        %----------------------------------------------------------
        % True nonlinear residual
        %----------------------------------------------------------
        fNew = assemble_RHS( ...
            P3, psiNew, quad, ...
            pprime, FFprime, dpsiNew);

        lhsI = ...
            K(freeNodes,:)*psiNew;

        residualVector = ...
            lhsI-fNew(freeNodes);

        residualScale = max([ ...
            norm(lhsI), ...
            norm(fNew(freeNodes)), ...
            eps]);

        relativeResidual = ...
            norm(residualVector)/residualScale;

        %----------------------------------------------------------
        % Normalized-flux diagnostic
        %
        % This expression is valid for general constant boundary
        % flux, although construct_RHS is currently restricted to
        % psiBoundary = 0.
        %----------------------------------------------------------
        psiNodal = ...
            (psiNew-psiAxisNew)/dpsiNew;

        minPsiN = min(psiNodal);
        maxPsiN = max(psiNodal);

        %----------------------------------------------------------
        % Store iteration history
        %----------------------------------------------------------
        updateHistory(iter) = ...
            relativeUpdate;

        residualHistory(iter) = ...
            relativeResidual;

        psiAxisHistory(iter) = ...
            psiAxisNew;

        dpsiHistory(iter) = ...
            dpsiNew;

        minPsiNHistory(iter) = ...
            minPsiN;

        maxPsiNHistory(iter) = ...
            maxPsiN;

        if opts.verbose
            fprintf( ...
                '%5d   %+14.7e   %+11.4e   %.4e   %.4e\n', ...
                iter, ...
                psiAxisNew, ...
                dpsiNew, ...
                relativeUpdate, ...
                relativeResidual);
        end

        %----------------------------------------------------------
        % Accept iteration
        %----------------------------------------------------------
        psi = psiNew;
        psiAxis = psiAxisNew;
        dpsi = dpsiNew;
        resolvedAxisMode = resolvedAxisModeNew;
        fCurrent = fNew;

        %----------------------------------------------------------
        % Convergence requires both criteria
        %----------------------------------------------------------
        if relativeUpdate < opts.updateTolerance && ...
                relativeResidual < opts.residualTolerance
            converged = true;
            break
        end
    end

    %==============================================================
    % Package result
    %==============================================================
    result = struct();

    result.converged = converged;
    result.iterations = iter;

    result.psiAxis = psiAxis;
    result.dpsi = dpsi;

    result.axisModeRequested = ...
        opts.axisMode;

    result.axisMode = ...
        resolvedAxisMode;

    result.boundaryValue = ...
        psiBoundary;

    result.initialBoundaryError = ...
        initialBoundaryError;

    result.symmetryError = ...
        symmetryError;

    result.operatorReused = true;

    result.updateHistory = ...
        updateHistory(1:iter);

    result.residualHistory = ...
        residualHistory(1:iter);

    result.psiAxisHistory = ...
        psiAxisHistory(1:iter);

    result.dpsiHistory = ...
        dpsiHistory(1:iter);

    result.minPsiNHistory = ...
        minPsiNHistory(1:iter);

    result.maxPsiNHistory = ...
        maxPsiNHistory(1:iter);

    result.finalRHS = ...
        fCurrent;

    result.finalResidual = ...
        residualVector;

    result.finalRelativeUpdate = ...
        updateHistory(iter);

    result.finalRelativeResidual = ...
        residualHistory(iter);

    result.freeNodes = ...
        freeNodes;

    result.boundaryNodes = ...
        boundaryNodes;

    result.profilePointCount = ...
        numel(pprime);

    %==============================================================
    % Report
    %==============================================================
    if opts.verbose
        fprintf('--------------------------------------------------------\n');

        if converged
            fprintf( ...
                'Picard iteration converged in %d iterations.\n', ...
                iter);
        else
            fprintf( ...
                ['Picard iteration did not converge in ', ...
                 '%d iterations.\n'], iter);
        end

        fprintf('Final psi_axis       = %+14.7e\n', ...
            psiAxis);

        fprintf('Final dpsi           = %+14.7e\n', ...
            dpsi);

        fprintf('Resolved axis mode   = %s\n', ...
            resolvedAxisMode);

        fprintf('Final update error   = %.4e\n', ...
            updateHistory(iter));

        fprintf('Final residual error = %.4e\n', ...
            residualHistory(iter));

        fprintf('Final nodal psiN     = [%.6e, %.6e]\n', ...
            minPsiNHistory(iter), ...
            maxPsiNHistory(iter));

        fprintf('========================================================\n\n');
    end

    if ~converged
        warning('GS:PicardPrepared:NoConvergence', ...
            ['Prepared Picard iteration reached the maximum ', ...
             'iteration count without satisfying both tolerances.']);
    end
end


function [psiAxis, dpsi, resolvedAxisMode] = ...
    evaluate_axis_and_span( ...
        psi, freeNodes, psiBoundary, ...
        axisMode, fluxSpanTolerance)
%EVALUATE_AXIS_AND_SPAN
%
% Determine the nodal magnetic-axis flux and signed flux span.

    freePsi = psi(freeNodes);

    switch axisMode
        case 'max'
            psiAxis = max(freePsi);
            resolvedAxisMode = 'max';

        case 'min'
            psiAxis = min(freePsi);
            resolvedAxisMode = 'min';

        case 'auto'
            psiMaximum = max(freePsi);
            psiMinimum = min(freePsi);

            maximumDistance = ...
                abs(psiMaximum-psiBoundary);

            minimumDistance = ...
                abs(psiMinimum-psiBoundary);

            if maximumDistance >= minimumDistance
                psiAxis = psiMaximum;
                resolvedAxisMode = 'max';
            else
                psiAxis = psiMinimum;
                resolvedAxisMode = 'min';
            end

        otherwise
            error('GS:PicardPrepared:AxisMode', ...
                'Unknown axis mode.');
    end

    dpsi = psiBoundary-psiAxis;

    fluxScale = max([ ...
        1, ...
        abs(psiBoundary), ...
        abs(psiAxis), ...
        norm(psi,inf)]);

    if abs(dpsi) <= ...
            fluxSpanTolerance*fluxScale
        error('GS:PicardPrepared:FluxSpan', ...
            ['psi_axis is too close to psi_boundary, ', ...
             'so normalized flux is undefined.']);
    end
end


function opts = parse_Picard_options(options)
%PARSE_PICARD_OPTIONS Validate prepared Picard options.

    if ~isstruct(options) || ~isscalar(options)
        error('GS:PicardPrepared:OptionsType', ...
            'options must be a scalar structure.');
    end

    defaults = struct();

    defaults.omega             = 0.5;
    defaults.maxIterations     = 100;
    defaults.updateTolerance   = 1.e-8;
    defaults.residualTolerance = 1.e-8;
    defaults.axisMode          = 'max';
    defaults.symmetryTolerance = 1.e-12;
    defaults.fluxSpanTolerance = 1.e-12;
    defaults.verbose           = true;

    unknownFields = ...
        setdiff(fieldnames(options),fieldnames(defaults));

    if ~isempty(unknownFields)
        error('GS:PicardPrepared:UnknownOption', ...
            'Unknown option: %s', ...
            strjoin(unknownFields,', '));
    end

    opts = defaults;

    optionNames = fieldnames(options);

    for k = 1:numel(optionNames)
        optionName = optionNames{k};
        opts.(optionName) = options.(optionName);
    end

    validateattributes(opts.omega, {'numeric'}, ...
        {'real', 'finite', 'scalar', '>', 0, '<=', 1});

    validateattributes(opts.maxIterations, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', 'positive'});

    validateattributes(opts.updateTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(opts.residualTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(opts.symmetryTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(opts.fluxSpanTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(opts.verbose, ...
        {'logical', 'numeric'}, ...
        {'real', 'finite', 'scalar'});

    opts.verbose = logical(opts.verbose);

    if isstring(opts.axisMode)
        if ~isscalar(opts.axisMode)
            error('GS:PicardPrepared:AxisMode', ...
                'axisMode must be scalar.');
        end

        opts.axisMode = char(opts.axisMode);
    end

    if ~ischar(opts.axisMode) || ...
            ~isrow(opts.axisMode)
        error('GS:PicardPrepared:AxisMode', ...
            'axisMode must be ''max'', ''min'', or ''auto''.');
    end

    opts.axisMode = lower(strtrim(opts.axisMode));

    if ~ismember(opts.axisMode, ...
            {'max', 'min', 'auto'})
        error('GS:PicardPrepared:AxisMode', ...
            'axisMode must be ''max'', ''min'', or ''auto''.');
    end
end


function validate_prepared_problem(problem)
%VALIDATE_PREPARED_PROBLEM Validate reusable FEM data.

    if ~isstruct(problem) || ~isscalar(problem)
        error('GS:PicardPrepared:ProblemType', ...
            'problem must be a scalar structure.');
    end

    requiredProblemFields = { ...
        'prepared', ...
        'mesh', ...
        'quad', ...
        'operator', ...
        'boundary'};

    assert_required_fields( ...
        problem, requiredProblemFields, 'problem');

    validateattributes(problem.prepared, ...
        {'logical', 'numeric'}, ...
        {'real', 'finite', 'scalar'});

    if ~logical(problem.prepared)
        error('GS:PicardPrepared:NotPrepared', ...
            'problem.prepared must be true.');
    end

    requiredOperatorFields = { ...
        'K', ...
        'Kfactor', ...
        'freeNodes', ...
        'boundaryNodes', ...
        'boundaryValues', ...
        'boundaryCorrection', ...
        'symmetryError'};

    assert_required_fields( ...
        problem.operator, ...
        requiredOperatorFields, ...
        'problem.operator');

    assert_required_fields( ...
        problem.boundary, ...
        {'psi'}, ...
        'problem.boundary');

    P3 = problem.mesh;
    Ndof = size(P3.points,1);

    K = problem.operator.K;

    if ~isnumeric(K) || ...
            ~isreal(K) || ...
            ~isequal(size(K),[Ndof,Ndof])
        error('GS:PicardPrepared:OperatorSize', ...
            'problem.operator.K must be a real Ndof-by-Ndof matrix.');
    end

    if any(~isfinite(nonzeros(K)))
        error('GS:PicardPrepared:InvalidOperator', ...
            'problem.operator.K contains NaN or Inf.');
    end

    if isempty(problem.operator.Kfactor)
        error('GS:PicardPrepared:MissingFactorization', ...
            'problem.operator.Kfactor is empty.');
    end

    freeNodes = problem.operator.freeNodes(:);
    boundaryNodes = problem.operator.boundaryNodes(:);

    if isempty(freeNodes) || isempty(boundaryNodes)
        error('GS:PicardPrepared:NodePartition', ...
            'Prepared free and boundary node sets must be nonempty.');
    end

    allNodes = sort([freeNodes; boundaryNodes]);

    if ~isequal(allNodes,(1:Ndof).')
        error('GS:PicardPrepared:NodePartition', ...
            ['Prepared free and boundary nodes do not form ', ...
             'a complete nonoverlapping partition.']);
    end

    if numel(unique(freeNodes)) ~= numel(freeNodes) || ...
            numel(unique(boundaryNodes)) ~= numel(boundaryNodes)
        error('GS:PicardPrepared:NodePartition', ...
            'Prepared node sets contain duplicated indices.');
    end

    if numel(problem.operator.boundaryValues) ~= ...
            numel(boundaryNodes)
        error('GS:PicardPrepared:BoundaryValues', ...
            'Prepared boundary values have incompatible size.');
    end

    if numel(problem.operator.boundaryCorrection) ~= ...
            numel(freeNodes)
        error('GS:PicardPrepared:BoundaryCorrection', ...
            'Prepared boundary correction has incompatible size.');
    end
end


function assert_required_fields( ...
    value, requiredFields, valueName)
%ASSERT_REQUIRED_FIELDS Check required structure fields.

    if ~isstruct(value) || ~isscalar(value)
        error('GS:PicardPrepared:StructureType', ...
            '%s must be a scalar structure.', valueName);
    end

    for k = 1:numel(requiredFields)
        fieldName = requiredFields{k};

        if ~isfield(value,fieldName)
            error('GS:PicardPrepared:MissingField', ...
                '%s.%s is required.', ...
                valueName,fieldName);
        end
    end
end