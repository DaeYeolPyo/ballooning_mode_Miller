function curves = salpha_marginal_contours(out)
%SALPHA_MARGINAL_CONTOURS Extract connected lambda_max=0 contour segments.
%
%   curves = salpha_marginal_contours(out)
%
% Returns a struct array with fields alpha, s, and n.  This uses MATLAB's
% contourc tracing on the full 2-D lambda grid, so it preserves connected
% marginal curves through folds better than row-by-row root extraction.

    alpha = out.alpha(:).';
    s = out.s(:).';
    lambdaPlot = out.lambda_max.';

    if size(lambdaPlot, 1) ~= numel(s) || size(lambdaPlot, 2) ~= numel(alpha)
        error('salpha_marginal_contours:BadSize', ...
            'lambda_max must have size numel(alpha)-by-numel(s).');
    end

    if min(lambdaPlot(:), [], 'omitnan') > 0 || max(lambdaPlot(:), [], 'omitnan') < 0
        curves = struct('alpha', {}, 's', {}, 'n', {});
        return;
    end

    C = contourc(alpha, s, lambdaPlot, [0 0]);
    curves = parse_contour_matrix(C);
end

function curves = parse_contour_matrix(C)
    curves = struct('alpha', {}, 's', {}, 'n', {});
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        cols = (k+1):(k+n);
        if cols(end) > size(C, 2)
            break;
        end

        curves(end+1).alpha = C(1, cols).'; %#ok<AGROW>
        curves(end).s = C(2, cols).';
        curves(end).n = n;
        k = k + n + 1;
    end
end
