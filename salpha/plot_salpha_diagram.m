function figures = plot_salpha_diagram(result, opts)
%PLOT_SALPHA_DIAGRAM Plot stability and maximizing ballooning phase.
%
%   figures = PLOT_SALPHA_DIAGRAM(result)
%
% result must be produced by scan_salpha_diagram.  The black lambda=0
% contour is the marginal-stability boundary; positive lambda is unstable.

    arguments
        result (1,1) struct
        opts.Visible (1,1) string = "on"
        opts.ShowPhase (1,1) logical = true
        opts.Title (1,1) string = ""
        opts.ContourLevels (1,1) double = 30
    end

    requiredFields = { ...
        'shat', 'alpha', 'theta0', 'lambdaMax', ...
        'bestTheta0', 'hasMarginalCrossing', ...
        'equilibriumShat', 'equilibriumAlpha', 'targetPsiN'};
    assert_fields(result, requiredFields);

    visible = lower(strtrim(opts.Visible));
    if ~ismember(visible, ["on", "off"])
        error('salpha:Plot:InvalidVisibility', ...
            'Visible must be "on" or "off".');
    end

    validateattributes(opts.ContourLevels, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2}, ...
        mfilename, 'ContourLevels');

    if strlength(opts.Title) == 0
        plotTitle = sprintf( ...
            'Direct flux-surface $s$-$\\alpha$, $\\psi_N=%.3f$', ...
            result.targetPsiN);
    else
        plotTitle = char(opts.Title);
    end

    stabilityFigure = figure( ...
        'Color', 'w', ...
        'Visible', char(visible), ...
        'Name', 's-alpha stability diagram');
    stabilityAxes = axes(stabilityFigure);

    contourf(stabilityAxes, ...
        result.alpha, result.shat, result.lambdaMax, ...
        opts.ContourLevels, 'LineStyle', 'none');
    hold(stabilityAxes, 'on');

    if result.hasMarginalCrossing
        contour(stabilityAxes, ...
            result.alpha, result.shat, result.lambdaMax, [0, 0], ...
            'k', 'LineWidth', 2.2);
    end

    equilibriumInRange = ...
        result.equilibriumAlpha >= min(result.alpha) && ...
        result.equilibriumAlpha <= max(result.alpha) && ...
        result.equilibriumShat >= min(result.shat) && ...
        result.equilibriumShat <= max(result.shat);

    if equilibriumInRange
        plot(stabilityAxes, ...
            result.equilibriumAlpha, result.equilibriumShat, ...
            'p', 'MarkerSize', 12, ...
            'MarkerFaceColor', [1.0, 0.85, 0.1], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.1);
    end

    colorbar(stabilityAxes);
    colormap(stabilityAxes, turbo);
    grid(stabilityAxes, 'on');
    box(stabilityAxes, 'on');
    stabilityAxes.Layer = 'top';
    stabilityAxes.TickDir = 'out';
    stabilityAxes.FontSize = 11;
    xlabel(stabilityAxes, '$\alpha$', 'Interpreter', 'latex');
    ylabel(stabilityAxes, '$\hat{s}$', 'Interpreter', 'latex');
    title(stabilityAxes, plotTitle, 'Interpreter', 'latex');

    figures = struct();
    figures.stability = stabilityFigure;
    figures.stabilityAxes = stabilityAxes;
    figures.phase = gobjects(0);
    figures.phaseAxes = gobjects(0);

    if ~opts.ShowPhase
        return;
    end

    phaseFigure = figure( ...
        'Color', 'w', ...
        'Visible', char(visible), ...
        'Name', 'Phase of maximum ballooning eigenvalue');
    phaseAxes = axes(phaseFigure);

    [alphaGrid, shatGrid] = meshgrid(result.alpha, result.shat);
    valid = isfinite(result.bestTheta0);

    scatter(phaseAxes, ...
        alphaGrid(valid), shatGrid(valid), 58, ...
        result.bestTheta0(valid)/pi, 's', 'filled', ...
        'MarkerEdgeColor', [0.18, 0.18, 0.18], ...
        'LineWidth', 0.25);

    colormap(phaseAxes, hsv(numel(result.theta0)));
    theta0OverPi = result.theta0/pi;
    if numel(theta0OverPi) > 1
        phaseHalfStep = pi/max(numel(theta0OverPi), 1)/pi;
    else
        phaseHalfStep = 1;
    end
    clim(phaseAxes, [ ...
        min(theta0OverPi)-phaseHalfStep, ...
        max(theta0OverPi)+phaseHalfStep]);

    phaseColorbar = colorbar(phaseAxes);
    phaseColorbar.Label.Interpreter = 'latex';
    phaseColorbar.Label.String = ...
        '$\theta_{0,\mathrm{max}}/\pi$';

    grid(phaseAxes, 'on');
    box(phaseAxes, 'on');
    phaseAxes.Layer = 'top';
    phaseAxes.TickDir = 'out';
    phaseAxes.FontSize = 11;
    xlabel(phaseAxes, '$\alpha$', 'Interpreter', 'latex');
    ylabel(phaseAxes, '$\hat{s}$', 'Interpreter', 'latex');
    title(phaseAxes, ...
        'Phase of maximum ballooning eigenvalue', ...
        'Interpreter', 'latex');

    figures.phase = phaseFigure;
    figures.phaseAxes = phaseAxes;
end


function assert_fields(value, names)
    if ~isstruct(value) || ~isscalar(value)
        error('salpha:Plot:InvalidResult', ...
            'result must be a scalar structure.');
    end

    for index = 1:numel(names)
        if ~isfield(value, names{index})
            error('salpha:Plot:MissingField', ...
                'result.%s is required.', names{index});
        end
    end
end
