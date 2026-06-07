function h = plot_salpha_diagram(out, varargin)
%PLOT_SALPHA_DIAGRAM Plot lambda_max and the marginal ideal-ballooning curve.

    ip = inputParser;
    addParameter(ip, 'Parent', [], @(x) isempty(x) || isgraphics(x));
    addParameter(ip, 'ShowLabels', true, @(x)islogical(x)&&isscalar(x));
    addParameter(ip, 'MarginalMode', 'contour', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'RootBranch', 'all', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'MinCurvePoints', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
    addParameter(ip, 'LineWidth', 2.0, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'Title', [], @(x) isempty(x) || ischar(x) || isstring(x));
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
    %lambdaPlot(lambdaPlot < 0.0) = 0.0;
    h = struct();
    h.contourf = contourf(ax, A, S, lambdaPlot, 32, 'LineStyle', 'none');
    hold(ax, 'on');
    h.colorbar = colorbar(ax);
    ylabel(h.colorbar, '$\lambda_{\max}$', 'Interpreter', 'latex');

    if min(lambdaPlot(:), [], 'omitnan') <= 0 && max(lambdaPlot(:), [], 'omitnan') >= 0
        mode = lower(string(opt.MarginalMode));
        if mode == "contour" || mode == "both"
            h.marginal = contour(ax, A, S, lambdaPlot, [0 0], ...
                'r--', 'LineWidth', opt.LineWidth);
        else
            h.marginal = [];
        end

        if mode == "curve" || mode == "curves" || mode == "both"
            curves = salpha_marginal_contours(out);
            h.marginal_curves = plot_contour_curves(ax, curves, ...
                opt.MinCurvePoints, opt.LineWidth);
        else
            h.marginal_curves = gobjects(0);
        end

        if mode == "roots"
            roots = salpha_marginal_roots(out);
            h.marginal_roots = plot_root_branches(ax, roots, ...
                lower(string(opt.RootBranch)), opt.LineWidth);
        else
            h.marginal_roots = gobjects(0);
        end

        validModes = ["contour", "curve", "curves", "roots", "both"];
        if ~any(mode == validModes)
            error('plot_salpha_diagram:BadMarginalMode', ...
                'Unknown MarginalMode: %s', mode);
        end
    else
        h.marginal = [];
        h.marginal_curves = gobjects(0);
        h.marginal_roots = gobjects(0);
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
        if ~isempty(opt.Title)
            titleText = opt.Title;
        elseif isfield(out, 'plot_title')
            titleText = out.plot_title;
        else
            titleText = 'Miller local equilibrium $s$-$\alpha$ diagram';
        end
        title(ax, titleText, 'Interpreter', 'latex');
    end
end

function hCurves = plot_contour_curves(ax, curves, minCurvePoints, lineWidth)
    keep = arrayfun(@(c) c.n >= minCurvePoints, curves);
    curves = curves(keep);
    hCurves = gobjects(1, numel(curves));

    for k = 1:numel(curves)
        hCurves(k) = plot(ax, curves(k).alpha, curves(k).s, 'r--', ...
            'LineWidth', lineWidth);
    end
end

function hRoots = plot_root_branches(ax, roots, branch, lineWidth)
    s = roots.s(:);
    branches = roots.alpha_branches;

    switch branch
        case "all"
            hRoots = gobjects(1, size(branches, 2));
            for k = 1:size(branches, 2)
                hRoots(k) = plot(ax, branches(:, k), s, 'r--', ...
                    'LineWidth', lineWidth);
            end

        case {"left", "leftmost"}
            hRoots = plot(ax, roots.alpha_left, s, 'r--', ...
                'LineWidth', lineWidth);

        case {"right", "rightmost"}
            hRoots = plot(ax, roots.alpha_right, s, 'r--', ...
                'LineWidth', lineWidth);

        otherwise
            error('plot_salpha_diagram:BadRootBranch', ...
                'Unknown RootBranch: %s', branch);
    end
end
