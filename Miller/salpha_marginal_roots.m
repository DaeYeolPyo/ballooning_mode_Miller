function roots = salpha_marginal_roots(out)
%SALPHA_MARGINAL_ROOTS Extract alpha roots of lambda_max(s, alpha)=0.
%
%   roots = salpha_marginal_roots(out)
%
% The input is the struct returned by scan_miller_salpha_diagram or the
% compatible C-shape scan. Roots are found by linear interpolation in alpha
% at each fixed s value.

    alpha = out.alpha(:);
    s = out.s(:);
    lambda = out.lambda_max;

    if size(lambda, 1) ~= numel(alpha) || size(lambda, 2) ~= numel(s)
        error('salpha_marginal_roots:BadSize', ...
            'lambda_max must have size numel(alpha)-by-numel(s).');
    end

    alphaByS = cell(numel(s), 1);
    nRoots = zeros(numel(s), 1);
    for is = 1:numel(s)
        alphaByS{is} = roots_for_profile(alpha, lambda(:, is));
        nRoots(is) = numel(alphaByS{is});
    end

    maxRoots = max(nRoots, [], 'omitnan');
    alphaBranches = nan(numel(s), maxRoots);
    for is = 1:numel(s)
        r = alphaByS{is};
        alphaBranches(is, 1:numel(r)) = r;
    end

    roots = struct();
    roots.s = s;
    roots.alpha_by_s = alphaByS;
    roots.alpha_branches = alphaBranches;
    roots.n_roots = nRoots;

    if maxRoots >= 1
        roots.alpha_left = alphaBranches(:, 1);
        roots.alpha_right = nan(numel(s), 1);
        for is = 1:numel(s)
            if nRoots(is) > 0
                roots.alpha_right(is) = alphaBranches(is, nRoots(is));
            end
        end
    else
        roots.alpha_left = nan(numel(s), 1);
        roots.alpha_right = nan(numel(s), 1);
    end
end

function roots = roots_for_profile(alpha, lambda)
    alpha = alpha(:);
    lambda = lambda(:);
    roots = [];

    for i = 1:numel(alpha)-1
        f1 = lambda(i);
        f2 = lambda(i+1);
        if ~isfinite(f1) || ~isfinite(f2)
            continue;
        end

        if f1 == 0
            roots(end+1, 1) = alpha(i); %#ok<AGROW>
        elseif f1*f2 < 0
            roots(end+1, 1) = alpha(i) ...
                - f1*(alpha(i+1)-alpha(i))/(f2-f1); %#ok<AGROW>
        end
    end

    if isfinite(lambda(end)) && lambda(end) == 0
        roots(end+1, 1) = alpha(end); %#ok<AGROW>
    end

    if ~isempty(roots)
        roots = unique(roots, 'stable');
    end
end
