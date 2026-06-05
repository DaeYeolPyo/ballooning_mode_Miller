function h = plot_miller_delta_marginal_curves(out, varargin)
%PLOT_MILLER_DELTA_MARGINAL_CURVES Plot lambda_max=0 curves for delta scans.

    ip = inputParser;
    addParameter(ip, 'Parent', [], @(x) isempty(x) || isgraphics(x));
    addParameter(ip, 'LineWidth', 2.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'ShowLegend', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    if isempty(opt.Parent)
        figure;
        ax = axes;
    else
        ax = opt.Parent;
    end
    hold(ax, 'on');

    nD = numel(out.delta);
    colors = lines(max(nD, 1));
    h = gobjects(1, nD);

    for id = 1:nD
        scan = out.scans{id};
        [A, S] = meshgrid(scan.alpha, scan.s);
        lam = scan.lambda_max.';

        hasCrossing = min(lam(:), [], 'omitnan') <= 0 && max(lam(:), [], 'omitnan') >= 0;
        if hasCrossing
            [~, h(id)] = contour(ax, A, S, lam, [0 0], ...
                'LineWidth', opt.LineWidth, ...
                'LineColor', colors(id,:), ...
                'DisplayName', sprintf('\\delta = %.3g', out.delta(id)));
        else
            warning('plot_miller_delta_marginal_curves:NoCrossing', ...
                'No lambda=0 crossing found for delta = %.6g.', out.delta(id));
        end
    end

    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$\alpha$', 'Interpreter', 'latex');
    ylabel(ax, '$\hat{s}$', 'Interpreter', 'latex');
    title(ax, 'Miller local equilibrium marginal stability curves', 'Interpreter', 'latex');

    if opt.ShowLegend
        legend(ax, 'Location', 'best');
    end
end
