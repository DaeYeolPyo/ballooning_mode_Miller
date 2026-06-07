function h = plot_miller_delta_marginal_curves(out, varargin)
%PLOT_MILLER_DELTA_MARGINAL_CURVES Plot lambda_max=0 curves for delta scans.

    ip = inputParser;
    addParameter(ip, 'Parent', [], @(x) isempty(x) || isgraphics(x));
    addParameter(ip, 'MarginalMode', 'contour', @(x)ischar(x)||isstring(x));
    addParameter(ip, 'MinCurvePoints', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
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
            mode = lower(string(opt.MarginalMode));
            switch mode
                case "contour"
                    [~, h(id)] = contour(ax, A, S, lam, [0 0], ...
                        'LineWidth', opt.LineWidth, ...
                        'LineColor', colors(id,:), ...
                        'DisplayName', sprintf('\\delta = %.3g', out.delta(id)));

                case {"curve", "curves"}
                    h(id) = plot_delta_curves(ax, scan, colors(id,:), ...
                        opt.MinCurvePoints, opt.LineWidth, ...
                        sprintf('\\delta = %.3g', out.delta(id)));

                otherwise
                    error('plot_miller_delta_marginal_curves:BadMarginalMode', ...
                        'Unknown MarginalMode: %s', mode);
            end
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

function hFirst = plot_delta_curves(ax, scan, color, minCurvePoints, lineWidth, displayName)
    curves = salpha_marginal_contours(scan);
    keep = arrayfun(@(c) c.n >= minCurvePoints, curves);
    curves = curves(keep);

    hFirst = gobjects(1);
    for k = 1:numel(curves)
        if k == 1
            name = displayName;
        else
            name = '';
        end

        hLine = plot(ax, curves(k).alpha, curves(k).s, ...
            'LineWidth', lineWidth, ...
            'Color', color, ...
            'DisplayName', name);

        if k == 1
            hFirst = hLine;
        end
    end
end
