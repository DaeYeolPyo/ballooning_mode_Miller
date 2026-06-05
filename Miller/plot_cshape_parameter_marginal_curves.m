function h = plot_cshape_parameter_marginal_curves(out, varargin)
%PLOT_CSHAPE_PARAMETER_MARGINAL_CURVES Plot lambda_max=0 curves.

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

    nP = numel(out.parameterValues);
    colors = lines(max(nP, 1));
    h = gobjects(1, nP);

    for iparam = 1:nP
        scan = out.scans{iparam};
        [A, S] = meshgrid(scan.alpha, scan.s);
        lam = scan.lambda_max.';

        hasCrossing = min(lam(:), [], 'omitnan') <= 0 ...
            && max(lam(:), [], 'omitnan') >= 0;
        if hasCrossing
            [~, h(iparam)] = contour(ax, A, S, lam, [0 0], ...
                'LineWidth', opt.LineWidth, ...
                'LineColor', colors(iparam,:), ...
                'DisplayName', sprintf('%s = %.4g', ...
                    out.parameterName, out.parameterValues(iparam)));
        else
            warning('plot_cshape_parameter_marginal_curves:NoCrossing', ...
                'No lambda=0 crossing found for %s = %.6g.', ...
                out.parameterName, out.parameterValues(iparam));
        end
    end

    grid(ax, 'on');
    box(ax, 'on');
    xlabel(ax, '$\alpha$', 'Interpreter', 'latex');
    ylabel(ax, '$\hat{s}$', 'Interpreter', 'latex');
    title(ax, sprintf('C-shape %s scan marginal stability curves', ...
        out.parameterName), 'Interpreter', 'latex');

    if opt.ShowLegend
        legend(ax, 'Location', 'best');
    end
end
