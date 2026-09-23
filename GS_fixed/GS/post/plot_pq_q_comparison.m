function figureHandle = plot_pq_q_comparison(equilibrium)
%PLOT_PQ_Q_COMPARISON Compare prescribed and final calculated q profiles.
%
%   figureHandle = plot_pq_q_comparison(equilibrium)
%
% This diagnostic is intended for profile mode 3 (prescribed p'-q).  It
% remains useful when the outer inverse iteration does not converge because
% solve_pq_constraint returns its last accepted, internally consistent
% equilibrium.

    arguments
        equilibrium (1,1) struct
    end

    assert_required_fields(equilibrium, ...
        {'qResult','profile','constraintResult'},'equilibrium');
    assert_required_fields(equilibrium.qResult, ...
        {'psiN','q'},'equilibrium.qResult');
    assert_required_fields(equilibrium.profile, ...
        {'qPrescribed'},'equilibrium.profile');
    assert_required_fields(equilibrium.profile.qPrescribed, ...
        {'psiN','q'},'equilibrium.profile.qPrescribed');

    calculatedPsiN = equilibrium.qResult.psiN(:);
    calculatedQ = equilibrium.qResult.q(:);
    prescribedPsiN = equilibrium.profile.qPrescribed.psiN(:);
    prescribedQ = equilibrium.profile.qPrescribed.q(:);

    validate_profile(calculatedPsiN,calculatedQ,'calculated');
    validate_profile(prescribedPsiN,prescribedQ,'prescribed');

    gridTolerance = 100*eps(max(1,max(abs( ...
        [calculatedPsiN;prescribedPsiN]))));

    if numel(calculatedPsiN) == numel(prescribedPsiN) && ...
            norm(calculatedPsiN-prescribedPsiN,inf) <= gridTolerance
        prescribedOnCalculatedGrid = prescribedQ;
    else
        if calculatedPsiN(1) < prescribedPsiN(1)-gridTolerance || ...
                calculatedPsiN(end) > prescribedPsiN(end)+gridTolerance
            error('GS:PlotPQ:QGridRange', ...
                ['The prescribed q grid does not cover the final ', ...
                 'calculated q grid.']);
        end

        prescribedOnCalculatedGrid = interp1( ...
            prescribedPsiN,prescribedQ,calculatedPsiN,'pchip');
    end

    if any(~isfinite(prescribedOnCalculatedGrid)) || ...
            any(prescribedOnCalculatedGrid == 0)
        error('GS:PlotPQ:InvalidInterpolatedTarget', ...
            'The prescribed q profile must be finite and nonzero.');
    end

    relativeError = ...
        (calculatedQ-prescribedOnCalculatedGrid) ...
        ./prescribedOnCalculatedGrid;
    [maximumError,maximumIndex] = max(abs(relativeError));

    result = equilibrium.constraintResult;
    converged = isfield(result,'converged') && ...
        isscalar(result.converged) && logical(result.converged);

    if converged
        statusText = 'converged';
    elseif isfield(result,'terminationReason') && ...
            ~isempty(result.terminationReason)
        statusText = sprintf( ...
            'not converged: %s',char(string(result.terminationReason)));
    else
        statusText = 'not converged';
    end

    figureHandle = figure( ...
        'Name','p-q constraint: target and calculated q', ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','normalized', ...
        'Position',[0.18,0.10,0.64,0.72]);

    layout = tiledlayout(figureHandle,2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    qAxes = nexttile(layout);
    plot(qAxes,prescribedPsiN,prescribedQ, ...
        'k--','LineWidth',1.8,'DisplayName','target q');
    hold(qAxes,'on');
    plot(qAxes,calculatedPsiN,calculatedQ, ...
        'Color',[0.85,0.15,0.10], ...
        'LineWidth',1.6, ...
        'DisplayName','final calculated q');
    grid(qAxes,'on');
    box(qAxes,'on');
    xlim(qAxes,[0,1]);
    ylabel(qAxes,'q');
    title(qAxes,'Prescribed and final safety-factor profiles');
    legend(qAxes,'Location','best');
    qAxes.FontSize = 11;

    errorAxes = nexttile(layout);
    plot(errorAxes,calculatedPsiN,100*relativeError, ...
        'b-','LineWidth',1.5, ...
        'DisplayName','(q_{calc}-q_{target})/q_{target}');
    hold(errorAxes,'on');
    yline(errorAxes,0,'k-','HandleVisibility','off');

    if isfield(result,'tolerances') && ...
            isstruct(result.tolerances) && ...
            isfield(result.tolerances,'qRelative') && ...
            isscalar(result.tolerances.qRelative) && ...
            isfinite(result.tolerances.qRelative) && ...
            result.tolerances.qRelative > 0
        tolerancePercent = 100*result.tolerances.qRelative;
        yline(errorAxes,tolerancePercent,'r:', ...
            'LineWidth',1.2, ...
            'DisplayName','maximum-error tolerance');
        yline(errorAxes,-tolerancePercent,'r:', ...
            'LineWidth',1.2, ...
            'HandleVisibility','off');
    end

    plot(errorAxes,calculatedPsiN(maximumIndex), ...
        100*relativeError(maximumIndex),'ko', ...
        'MarkerFaceColor','y', ...
        'DisplayName','largest plotted error');
    grid(errorAxes,'on');
    box(errorAxes,'on');
    xlim(errorAxes,[0,1]);
    xlabel(errorAxes,'\psi_N');
    ylabel(errorAxes,'Relative q error [%]');
    title(errorAxes,sprintf( ...
        'Largest plotted error = %.3f%% at \\psi_N = %.3f', ...
        100*maximumError,calculatedPsiN(maximumIndex)));
    legend(errorAxes,'Location','best');
    errorAxes.FontSize = 11;

    title(layout,sprintf('p''-q constraint diagnostic (%s)',statusText), ...
        'Interpreter','none');
end


function validate_profile(psiN,q,profileName)
    if numel(psiN) ~= numel(q) || numel(psiN) < 2
        error('GS:PlotPQ:ProfileSize', ...
            'The %s q profile has an invalid size.',profileName);
    end

    if any(~isfinite(psiN)) || any(diff(psiN) <= 0) || ...
            any(~isfinite(q))
        error('GS:PlotPQ:InvalidProfile', ...
            ['The %s q profile must have a finite increasing psiN ', ...
             'grid and finite q values.'],profileName);
    end
end


function assert_required_fields(value,requiredFields,valueName)
    if ~isstruct(value) || ~isscalar(value)
        error('GS:PlotPQ:InvalidStructure', ...
            '%s must be a scalar structure.',valueName);
    end

    for k = 1:numel(requiredFields)
        if ~isfield(value,requiredFields{k})
            error('GS:PlotPQ:MissingField', ...
                '%s.%s is required.',valueName,requiredFields{k});
        end
    end
end
