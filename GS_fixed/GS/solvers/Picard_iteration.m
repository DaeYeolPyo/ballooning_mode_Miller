function [psi, result] = Picard_iteration(P3, quad, pprime, FFprime, psi0, options)
% [PICARD_ITERATION]
%
% Solve the nonlinear fixed-boundary Grad-Shafranov equation
% using under-relaxed Picard iteration.
%
% Current boundary convention:
%
%       psi_boundary = 0
%
% Current profile convention:
%
%       dp_dpsiN  = dp/dpsi_N
%       FdF_dpsiN = F*dF/dpsi_N
%
% where
%
%       psi_N = (psi-psi_axis)/(psi_boundary-psi_axis).
%
%
% INPUT
%
%   P3
%       P3 mesh structure.
%
%   quad
%       Precomputed P3 quadrature structure.
%
%   dp_dpsiN
%       Profile values of dp/dpsi_N on a uniform psi_N grid.
%
%   FdF_dpsiN
%       Profile values of F*dF/dpsi_N on the same grid.
%
%   psi0
%       Initial dimensional poloidal-flux vector.
%
%   options
%       Optional scalar structure:
%
%           omega               Picard relaxation parameter
%                               default: 0.5
%
%           maxIterations       maximum number of iterations
%                               default: 100
%
%           updateTolerance     relative update tolerance
%                               default: 1e-8
%
%           residualTolerance   relative nonlinear residual
%                               default: 1e-8
%
%           axisMode            'max', 'min', or 'auto'
%                               default: 'max'
%
%           symmetryTolerance   accepted KII symmetry error
%                               default: 1e-12
%
%           fluxSpanTolerance   relative minimum |dpsi|
%                               default: 1e-12
%
%           verbose             print iteration history
%                               default: true
%
%
% OUTPUT
%
%   psi
%       Converged dimensional poloidal-flux vector.
%
%   result
%       Convergence information and iteration histories.
    narginchk(5, 6);
    if nargin < 6 || isempty(options)
        options = struct();
    end

    opts = parse_Picard_options(options);

    %======================================================
    % Input validation
    %======================================================
    Ndof = size(P3.points, 1);
    psi = psi0(:);
    psiBoundary = 0.0;

    %======================================================
    % Assemble K matrix
    %======================================================
    K = assemble_K_matrix(P3, quad);

    if ~issparse(K)
        warning('[PICARD_ITERATION] K is not sparse.');
    end

    % Obtain boudnary node sets
    [KII, ~, bc] = apply_BC(P3, K, zeros(Ndof, 1), psiBoundary);

    freeNodes     = bc.freeNodes;
    boundaryNodes = bc.boundaryNodes;
    psiB          = bc.boundaryValues;

    % Enforce the boundary value exactly in the initial guess
    initialBoundaryError = norm(psi(boundaryNodes) - psiB, inf);

    psi(boundaryNodes) = psiB;

    %======================================================
    % Check symmetry before using one triangle of KII
    %======================================================
    symmetryError = norm(KII - KII.', 'fro')/max(norm(KII, 'fro'), eps);

    if symmetryError > opts.symmetryTolerance
        error('[PICARD_ITERATION] KII symmetry error %.3e exceeds tolerance %.3e', ...
            symmetryError, opts.symmetryTolerance);
    end

    % KIB and its fixed bounary correction do not change.
    KIB = K(freeNodes, boundaryNodes);
    boundaryCorrection = KIB*psiB;

    %======================================================
    % Factorize once and reuse during every Picard iteration
    %
    % 'lower' is used because KII may differ from KII' by
    % roundoff-sized amounts even after passing the symmetry test.
    %======================================================
    try
        Kfactor = decomposition(KII, 'chol', 'lower');
    catch ME
        error(['[PICARD_ITERATION] Failed to construct the Cholesky' ...
            'factorization of KII:\n%s'], ME.message);
    end

    %======================================================
    % Initial source evaluation
    %======================================================
    [psiAxis, dpsi] = evaluate_axis_and_span(psi, freeNodes, ...
        psiBoundary, opts.axisMode, opts.fluxSpanTolerance);

    fCurrent = assemble_RHS(P3, psi, quad, pprime, FFprime, dpsi);

    %======================================================
    % History arrays
    %======================================================
    maxIter = opts.maxIterations;

    updateHistory   = nan(maxIter, 1);
    residualHistory = nan(maxIter, 1);
    psiAxisHistory  = nan(maxIter, 1);
    dpsiHistory     = nan(maxIter, 1);
    minPsiNHistory   = nan(maxIter, 1);
    maxPsiNHistory   = nan(maxIter, 1);

    converged = false;

    if opts.verbose
        fprintf('\n');
        fprintf('=======================================================\n');
        fprintf(' Grad-Shafranov Picard iteration\n');
        fprintf('=======================================================\n');
        fprintf(' iter     psi_axis     dpsi     update error     residual error\n');
    end

    %======================================================
    % Picard iteration main loop
    %======================================================
    for iter = 1:maxIter
        %--------------------------------------------------
        % Solve using source evaluated at the current iteration
        %--------------------------------------------------
        rhsI = fCurrent(freeNodes) - boundaryCorrection;

        psiSolved = psi;

        psiSolved(freeNodes) = Kfactor\rhsI;
        psiSolved(boundaryNodes) = psiB;

        %--------------------------------------------------
        % Under-relaxation on free DOFs only
        %--------------------------------------------------
        psiNew = psi;

        psiNew(freeNodes) = (1 - opts.omega)*psi(freeNodes) ...
            + opts.omega*psiSolved(freeNodes);

        % Prevent relaxation or roundoff from changing the LCFS
        psiNew(boundaryNodes) = psiB;

        %--------------------------------------------------
        % Relative solution-update norm
        %--------------------------------------------------
        updateVector = psiNew(freeNodes) - psi(freeNodes);

        updateScale = max([norm(psiNew(freeNodes)), ...
            norm(psi(freeNodes)), eps]);

        relativeUpdate = norm(updateVector)/updateScale;

        %--------------------------------------------------
        % Recalculate axis and dpsi at the updated iteration
        %--------------------------------------------------
        [psiAxisNew, dpsiNew] = evaluate_axis_and_span(...
            psiNew, freeNodes, psiBoundary, ...
            opts.axisMode, opts.fluxSpanTolerance);

        %--------------------------------------------------
        % Evaluate the true nonlinear residual at psiNew
        %
        % This is not the lagged Picard residual.
        % The RHS is evaluated again using psiNew and dpsiNew.
        %--------------------------------------------------
        fNew = assemble_RHS(P3, psiNew, quad, pprime, FFprime, dpsiNew);

        lhsI = K(freeNodes, :)*psiNew;
        residualVector = lhsI - fNew(freeNodes);

        residualScale = max([norm(lhsI), norm(fNew(freeNodes)), eps]);

        relativeResidual = norm(residualVector)/residualScale;

        %--------------------------------------------------
        % Nodal normalized-flux range diagnostic
        %
        % For psi_boundary = 0:
        %   psi_N = (psi + dpsi)/dpsi
        %--------------------------------------------------
        psiNodal = (psiNew + dpsiNew)/dpsiNew;

        minPsiN = min(psiNodal);
        maxPsiN = max(psiNodal);

        %--------------------------------------------------
        % Store history
        %--------------------------------------------------
        updateHistory(iter)   = relativeUpdate;
        residualHistory(iter) = relativeResidual;
        psiAxisHistory(iter)  = psiAxisNew;
        dpsiHistory(iter)     = dpsiNew;
        minPsiNHistory(iter)  = minPsiN;
        maxPsiNHistory(iter)  = maxPsiN;

        if opts.verbose
            fprintf('%5d   %+14.7e   %+11.4e   %.4e   %.4e\n', ...
                iter, psiAxisNew, dpsiNew, ...
                relativeUpdate, relativeResidual);
        end

        % Accept the relaxded iteration
        psi      = psiNew;
        psiAxis  = psiAxisNew;
        dpsi     = dpsiNew;
        fCurrent = fNew;

        %--------------------------------------------------
        % Convergence requires both conditions
        %--------------------------------------------------
        if relativeUpdate < opts.updateTolerance && ...
            relativeResidual < opts.residualTolerance
            converged = true;
            break
        end
    end

    %======================================================
    % Packge result
    %======================================================
    result = struct();

    result.converged            = converged;
    result.iterations           = iter;
    result.psiAxis              = psiAxis;
    result.dpsi                 = dpsi;
    result.boundaryValue        = psiBoundary;
    result.initialBoundaryError = initialBoundaryError;
    result.symmetryError        = symmetryError;

    result.updateHistory   = updateHistory(1:iter);
    result.residualHistory = residualHistory(1:iter);
    result.psiAxisHistory  = psiAxisHistory(1:iter);
    result.dpsiHistory     = dpsiHistory(1:iter);
    result.minPsiNHistory  = minPsiNHistory(1:iter);
    result.maxPsiNHistory  = maxPsiNHistory(1:iter);

    result.finalRHS      = fCurrent;
    result.finalResidual = residualVector;
    result.freeNodes     = freeNodes;
    result.boundaryNodes = boundaryNodes;

    if opts.verbose
        fprintf('--------------------------------------------------------\n');

        if converged
            fprintf('Picard iteration converged in %d iterations.\n',iter);
        else
            fprintf('Picard iteration did not converge in %d iterations.\n', ...
                iter);
        end

        fprintf('Final psi_axis       = %+14.7e\n',psiAxis);
        fprintf('Final dpsi           = %+14.7e\n',dpsi);
        fprintf('Final update error   = %.4e\n', ...
            updateHistory(iter));
        fprintf('Final residual error = %.4e\n', ...
            residualHistory(iter));
        fprintf('Final nodal psiN     = [%.6e, %.6e]\n', ...
            minPsiNHistory(iter),maxPsiNHistory(iter));
        fprintf('========================================================\n\n');
    end

    if ~converged
        warning('GS:Picard:NoConvergence', ...
            ['Picard iteration reached the maximum iteration ', ...
             'count without satisfying both tolerances.']);
    end
end

function [psiAxis, dpsi] = evaluate_axis_and_span( ...
    psi, freeNodes, psiBoundary, axisMode, fluxSpanTolerance)
    % Determine the nodal magnetic-axis flux and signed flux span
    freePsi = psi(freeNodes);

    switch axisMode
        case 'max'
            psiAxis = max(freePsi);
        case 'min'
            psiAxis = min(freePsi);
        case 'auto'
            psiMaximum = max(freePsi);
            psiMinimum = min(freePsi);

            maximumDistance = ...
                abs(psiMaximum-psiBoundary);

            minimumDistance = ...
                abs(psiMinimum-psiBoundary);

            if maximumDistance >= minimumDistance
                psiAxis = psiMaximum;
            else
                psiAxis = psiMinimum;
            end

        otherwise
            error('[PICARD_ITERATION] Unknown axis mode.');
    end

    dpsi = psiBoundary - psiAxis;

    fluxScale = max([1, abs(psiBoundary), abs(psiAxis), norm(psi, inf)]);

    if abs(dpsi) <= fluxSpanTolerance*fluxScale
        error(['[PICARD_ITERATION] psi_axis is too close to', ...
            'psi_boundary, so normalized flux is undefined.']);
    end
end


function opts = parse_Picard_options(options)
    if ~isstruct(options) || ~isscalar(options)
        error('[PICARD_ITERATION] options must be a scalar structure.');
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

    unknownFields = setdiff(fieldnames(options), fieldnames(defaults));

    if ~isempty(unknownFields)
        error('[PICARD_ITERATION] Unknown optino: %s', ...
            strjoin(unknownFields, ', '));
    end

    opts = defaults;

    optionNames = fieldnames(options);

    for k = 1:numel(optionNames)
        name = optionNames{k};
        opts.(name) = options.(name);
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
    validateattributes(opts.verbose, {'logical', 'numeric'}, ...
        {'real', 'finite', 'scalar'});

    opts.verbose = logical(opts.verbose);

    if isstring(opts.axisMode)
        if ~isscalar(opts.axisMode)
            error('[PICARD_ITERATION] axisMode must be scalar.');
        end

        opts.axisMode = char(opts.axisMode);
    end

    if ~ischar(opts.axisMode) || ~isrow(opts.axisMode)
        error('[PICARD_ITERATION] axisMode must be ''max'', ''min'', or ''auto''.');
    end

    opts.axisMode = lower(strtrim(opts.axisMode));

    if ~ismember(opts.axisMode, {'max', 'min', 'auto'})
        error('[PICARD_ITERATION] axisMode must be ''max'', ''min'', or ''auto''.');
    end
end