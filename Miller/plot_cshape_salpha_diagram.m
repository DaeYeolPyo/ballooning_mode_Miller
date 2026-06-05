function h = plot_cshape_salpha_diagram(out, varargin)
%PLOT_CSHAPE_SALPHA_DIAGRAM Plot C-shape lambda_max and marginal curve.

    ip = inputParser;
    addParameter(ip, 'Parent', [], @(x) isempty(x) || isgraphics(x));
    addParameter(ip, 'ShowLabels', true, @(x)islogical(x)&&isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    if isempty(opt.Parent)
        figure;
        ax = axes;
    else
        ax = opt.Parent;
    end

    axes(ax);
    [A, S] = meshgrid(out.alpha, out.s);
    lambdaPlot = out.lambda_max.';

    h = struct();
    h.contourf = contourf(ax, A, S, lambdaPlot, 32, 'LineStyle', 'none');
    hold(ax, 'on');
    h.colorbar = colorbar(ax);
    ylabel(h.colorbar, '$\lambda_{\max}$', 'Interpreter', 'latex');

    if min(lambdaPlot(:), [], 'omitnan') <= 0 && max(lambdaPlot(:), [], 'omitnan') >= 0
        h.marginal = contour(ax, A, S, lambdaPlot, [0 0], ...
            'k--', 'LineWidth', 2.0);
    else
        h.marginal = [];
    end

    try
        colormap(ax, turbo);
    catch
        colormap(ax, parula);
    end
    grid(ax, 'on');
    box(ax, 'on');

    if opt.ShowLabels
        xlabel(ax, '$\alpha$', 'Interpreter', 'latex');
        ylabel(ax, '$\hat{s}$', 'Interpreter', 'latex');
        title(ax, sprintf('C-shape local equilibrium $s$-$\\alpha$ diagram, $\\psi_N=%.3g$', ...
            out.psin), 'Interpreter', 'latex');
    end
end
