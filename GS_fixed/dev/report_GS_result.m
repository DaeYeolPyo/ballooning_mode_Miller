function report_GS_result(equilibrium)
%[REPORT_GS_RESULT]
% Print a summary of a completed Grad-Shafranov solve.
%
%   report_GS_result(equilibrium)
%
% The input must be the structure returned by solve_GS.  This routine
% reports only quantities available from the numerical solution; it does
% not print verification errors that require an exact equilibrium.

    arguments
        equilibrium (1,1) struct
    end

    requiredEquilibriumFields = { ...
        'mesh', 'psi', 'psiN', 'axisNode', 'axisPoint', 'solverResult'};

    assert_required_fields( ...
        equilibrium, requiredEquilibriumFields, 'equilibrium');

    result = equilibrium.solverResult;

    requiredResultFields = { ...
        'converged', 'iterations', 'psiAxis', 'dpsi', ...
        'boundaryValue', 'updateHistory', 'residualHistory', ...
        'boundaryNodes', 'symmetryError'};

    assert_required_fields(result, requiredResultFields, ...
        'equilibrium.solverResult');

    P3 = equilibrium.mesh;

    requiredMeshFields = { ...
        'points', 'elements', 'P1points', 'boundaryNodes'};

    assert_required_fields(P3, requiredMeshFields, 'equilibrium.mesh');

    psi = equilibrium.psi(:);
    psiN = equilibrium.psiN(:);
    boundaryNodes = result.boundaryNodes(:);

    if isempty(result.updateHistory) || isempty(result.residualHistory)
        error(['[REPORT_GS_RESULT] The Picard histories must contain ', ...
            'at least one iteration.']);
    end

    if numel(psi) ~= size(P3.points, 1)
        error(['[REPORT_GS_RESULT] equilibrium.psi must contain one ', ...
            'value per P3 node.']);
    end

    if numel(psiN) ~= numel(psi)
        error(['[REPORT_GS_RESULT] equilibrium.psiN and equilibrium.psi ', ...
            'must have equal lengths.']);
    end

    if any(~isfinite(psi)) || any(~isfinite(psiN))
        error('[REPORT_GS_RESULT] The reported flux contains NaN or Inf.');
    end

    if isempty(boundaryNodes) || any(boundaryNodes < 1) || ...
            any(boundaryNodes > numel(psi)) || ...
            any(boundaryNodes ~= round(boundaryNodes))
        error('[REPORT_GS_RESULT] Invalid boundary-node indices.');
    end

    boundaryValues = expand_boundary_values( ...
        result.boundaryValue, boundaryNodes, numel(psi));

    finalBoundaryError = norm( ...
        psi(boundaryNodes)-boundaryValues, inf);

    finalUpdate = result.updateHistory(end);
    finalResidual = result.residualHistory(end);

    psiNMinimum = min(psiN);
    psiNMaximum = max(psiN);

    axisPoint = equilibrium.axisPoint(:).';

    if numel(axisPoint) ~= 2 || any(~isfinite(axisPoint))
        error(['[REPORT_GS_RESULT] equilibrium.axisPoint must contain ', ...
            'finite [R,Z] coordinates.']);
    end

    if result.psiAxis > mean(boundaryValues)
        axisExtremum = 'max';
    else
        axisExtremum = 'min';
    end

    if logical(result.converged)
        solveStatus = 'CONVERGED';
    else
        solveStatus = 'NOT CONVERGED';
    end

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf(' Fixed-boundary Grad-Shafranov solve report\n');
    fprintf('============================================================\n');
    fprintf('Status                     = %s\n', solveStatus);
    fprintf('Picard iterations          = %d\n', result.iterations);
    fprintf('Final relative update      = %.3e\n', finalUpdate);
    fprintf('Final relative residual    = %.3e\n', finalResidual);
    fprintf('\n');
    fprintf('psi_axis                   = %+14.7e\n', result.psiAxis);
    fprintf('psi_boundary               = %+14.7e\n', ...
        mean(boundaryValues));
    fprintf('dpsi                       = %+14.7e\n', result.dpsi);
    fprintf('Nodal psi_N range          = [%+.6e, %+.6e]\n', ...
        psiNMinimum, psiNMaximum);
    fprintf('Axis extremum              = %s\n', axisExtremum);
    fprintf('Nodal magnetic axis        = node %d at (%.8f, %.8f)\n', ...
        equilibrium.axisNode, axisPoint(1), axisPoint(2));
    fprintf('Final boundary error       = %.3e\n', finalBoundaryError);
    fprintf('Reduced-matrix symmetry    = %.3e\n', result.symmetryError);

    fprintf('\n');
    fprintf('P1 vertices                = %d\n', size(P3.P1points, 1));
    fprintf('P1 triangles               = %d\n', size(P3.elements, 1));
    fprintf('P3 DOFs                    = %d\n', size(P3.points, 1));
    fprintf('P3 boundary DOFs           = %d\n', numel(P3.boundaryNodes));
    
    fprintf('\n');
    fprintf('Plasma current            = %.6e MA\n', ...
        equilibrium.global.Ip*1.e-6);
    fprintf('Plasma volume             = %.6e m^3\n', ...
        equilibrium.global.volume);
    fprintf('Toroidal beta             = %.6e  (%.4f %%)\n', ...
        equilibrium.global.betaTor, ...
        100*equilibrium.global.betaTor);
    fprintf('Poloidal beta             = %.6e\n', ...
        equilibrium.global.betaPol);
    fprintf('Normalized beta           = %.6e\n', ...
        equilibrium.global.betaN);
    fprintf('Internal inductance li_3  = %.6e\n', ...
        equilibrium.global.li3);
    fprintf('============================================================\n');
end


function assert_required_fields(value, requiredFields, valueName)
    missingFields = setdiff(requiredFields, fieldnames(value));

    if ~isempty(missingFields)
        error('[REPORT_GS_RESULT] %s is missing field(s): %s.', ...
            valueName, strjoin(missingFields, ', '));
    end
end


function boundaryValues = expand_boundary_values( ...
    boundaryValue, boundaryNodes, Ndof)

    if ~isnumeric(boundaryValue) || ~isreal(boundaryValue) || ...
            any(~isfinite(boundaryValue(:)))
        error(['[REPORT_GS_RESULT] solverResult.boundaryValue must be ', ...
            'finite and real.']);
    end

    Nb = numel(boundaryNodes);

    if isscalar(boundaryValue)
        boundaryValues = repmat(boundaryValue, Nb, 1);
    elseif isvector(boundaryValue) && numel(boundaryValue) == Nb
        boundaryValues = boundaryValue(:);
    elseif isvector(boundaryValue) && numel(boundaryValue) == Ndof
        fullBoundaryValues = boundaryValue(:);
        boundaryValues = fullBoundaryValues(boundaryNodes);
    else
        error(['[REPORT_GS_RESULT] solverResult.boundaryValue must be a ', ...
            'scalar, a boundary-node vector, or a full nodal vector.']);
    end
end
