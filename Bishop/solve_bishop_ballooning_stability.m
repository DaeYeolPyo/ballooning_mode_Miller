function stab = solve_bishop_ballooning_stability(ch3, varargin)
%SOLVE_BISHOP_BALLOONING_STABILITY Solve the Bishop ballooning stability problem.
%
%   stab = SOLVE_BISHOP_BALLOONING_STABILITY(ch3)
%   stab = SOLVE_BISHOP_BALLOONING_STABILITY(ch3, 'Name', value, ...)
%
% Input
%   ch3 : output of compute_bishop_ch3_terms(...). If ch3.extended exists,
%         the extended ballooning-theta domain is used by default.
%
% Name-value options
%   'UseExtended'       : true/false, default true
%   'ThetaMax'          : optional finite theta cutoff. If supplied, only
%                         points with abs(theta) <= ThetaMax are retained.
%   'NumModes'          : number of eigenmodes. Default 6.
%   'UseEigs'           : true/false, default true.
%   'Target'            : eigs target. Default 'auto'. In auto mode the
%                         lowest eigenvalue is first bracketed by Cholesky
%                         tests, then eigs uses a nearby numeric shift.
%   'StabilityTolerance': eigenvalue tolerance around zero. Default 1e-8.
%   'Plot'              : true/false. Default false.
%
% The operator is assembled in the divided, self-adjoint form
%
%   H F = lambda F,
%   H = - d/dl ( A_inner d/dl ) - V,
%   V = C_drive / Bp.
%
% The corresponding energy functional is
%
%   deltaW[F] = int ( A_inner*(dF/dl)^2 - V*F^2 ) dl.
%
% Stability convention
%   lambda_min > 0 : stable
%   lambda_min = 0 : marginal
%   lambda_min < 0 : unstable
%
% Boundary condition
%   F = 0 at the two ends of the finite extended-theta window.

    p = inputParser;
    addParameter(p, 'UseExtended', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'ThetaMax', [], @(x)isnumeric(x) && isscalar(x) && x > 0 || isempty(x));
    addParameter(p, 'NumModes', 6, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'UseEigs', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Target', 'auto', @(x)ischar(x) || isstring(x));
    addParameter(p, 'StabilityTolerance', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Plot', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    path = choose_path(ch3, opt);
    path = restrict_theta_window(path, opt.ThetaMax);

    mats = assemble_energy_matrices(path);
    lambdaMin = estimate_lowest_eigenvalue(mats.H, mats.M);
    sol = solve_energy_eigs(mats, opt, lambdaMin);

    tol = opt.StabilityTolerance;
    if lambdaMin < -tol
        status = 'unstable';
    elseif lambdaMin > tol
        status = 'stable';
    else
        status = 'marginal';
    end

    energies = compute_mode_energies(mats, sol);

    stab = struct();
    stab.status = status;
    stab.lambdaMin = lambdaMin;
    stab.lambda = sol.lambda;
    stab.theta = mats.theta;
    stab.l = mats.l;
    stab.F = sol.F;
    stab.Fi = sol.Fi;
    stab.energy = energies;
    stab.mats = mats;
    stab.options = opt;

    if opt.Plot
        plot_stability_solution(stab);
    end
end

%======================================================================
function path = choose_path(ch3, opt)
    if opt.UseExtended && isfield(ch3, 'extended') && ~isempty(ch3.extended)
        path = ch3.extended;
        path.domain = 'extended';
    else
        path = ch3;
        path.domain = 'one-period';
    end

    if ~isfield(path, 'eq32') || ~isfield(path.eq32, 'A_inner')
        error('solve_bishop_ballooning_stability:MissingEq32', ...
            'Input must contain eq32.A_inner from compute_bishop_ch3_terms.');
    end

    if ~isfield(path.eq32, 'C_dividedByBp')
        if isfield(path.eq32, 'C_drive') && isfield(path.eq32, 'outerBp')
            path.eq32.C_dividedByBp = path.eq32.C_drive ./ path.eq32.outerBp;
        else
            error('solve_bishop_ballooning_stability:MissingDrive', ...
                'Input must contain eq32.C_dividedByBp or C_drive/outerBp.');
        end
    end
end

%======================================================================
function path = restrict_theta_window(path, thetaMax)
    if isempty(thetaMax)
        return;
    end

    theta = path.theta(:);
    keep = abs(theta) <= thetaMax;

    if nnz(keep) < 5
        error('solve_bishop_ballooning_stability:ThetaWindowTooSmall', ...
            'ThetaMax leaves fewer than five grid points.');
    end

    path.theta = path.theta(keep);
    path.l = path.l(keep);
    path.eq32.A_inner = path.eq32.A_inner(keep);
    path.eq32.C_dividedByBp = path.eq32.C_dividedByBp(keep);

    if isfield(path.eq32, 'C_drive')
        path.eq32.C_drive = path.eq32.C_drive(keep);
    end
    if isfield(path.eq32, 'outerBp')
        path.eq32.outerBp = path.eq32.outerBp(keep);
    end
end

%======================================================================
function mats = assemble_energy_matrices(path)
    l = path.l(:);
    theta = path.theta(:);
    Acoef = path.eq32.A_inner(:);
    Vcoef = path.eq32.C_dividedByBp(:);

    n = numel(l);
    if n < 5
        error('solve_bishop_ballooning_stability:TooFewPoints', ...
            'Need at least five points for Dirichlet stability solve.');
    end

    if any(diff(l) <= 0)
        error('solve_bishop_ballooning_stability:BadGrid', ...
            'The l grid must be strictly increasing.');
    end

    if any(~isfinite(Acoef)) || any(~isfinite(Vcoef))
        error('solve_bishop_ballooning_stability:NonFiniteCoefficients', ...
            'A_inner and C_dividedByBp must be finite.');
    end

    if any(Acoef <= 0)
        warning('solve_bishop_ballooning_stability:NonPositiveA', ...
            'A_inner has non-positive values; the energy operator may be ill-posed.');
    end

    Kfull = spalloc(n, n, 3*n);
    Pfull = spalloc(n, n, 3*n);
    Mfull = spalloc(n, n, 3*n);

    for e = 1:n-1
        h = l(e+1) - l(e);
        Ae = 0.5 * (Acoef(e) + Acoef(e+1));
        Ve = 0.5 * (Vcoef(e) + Vcoef(e+1));

        Ke = (Ae / h) .* [1 -1; -1 1];
        Me = (h / 6) .* [2 1; 1 2];
        Pe = Ve .* Me;

        nodes = [e, e+1];
        Kfull(nodes, nodes) = Kfull(nodes, nodes) + Ke;
        Mfull(nodes, nodes) = Mfull(nodes, nodes) + Me;
        Pfull(nodes, nodes) = Pfull(nodes, nodes) + Pe;
    end

    Hfull = Kfull - Pfull;
    interior = 2:n-1;

    mats = struct();
    mats.l = l;
    mats.theta = theta;
    mats.Acoef = Acoef;
    mats.Vcoef = Vcoef;
    mats.Kfull = Kfull;
    mats.Pfull = Pfull;
    mats.Mfull = Mfull;
    mats.Hfull = Hfull;
    mats.K = Kfull(interior, interior);
    mats.P = Pfull(interior, interior);
    mats.M = Mfull(interior, interior);
    mats.H = Hfull(interior, interior);
    mats.interior = interior;
    mats.n = n;
    mats.ni = numel(interior);
    mats.domain = path.domain;
end

%======================================================================
function sol = solve_energy_eigs(mats, opt, lambdaMinEstimate)
    nmodes = min(round(opt.NumModes), mats.ni);

    useEigs = opt.UseEigs && mats.ni > nmodes + 2;

    if useEigs
        try
            target = choose_eigs_target(opt.Target, lambdaMinEstimate);
            [Fi, D] = eigs(mats.H, mats.M, nmodes, target);
            lambda = diag(D);
            if any(~isfinite(lambda))
                error('solve_bishop_ballooning_stability:NonFiniteEigs', ...
                    'eigs returned non-finite eigenvalues.');
            end
        catch ME
            warning('solve_bishop_ballooning_stability:EigsFailed', ...
                'eigs failed (%s). Trying smallestabs/full fallback.', ME.message);
            try
                [Fi, D] = eigs(mats.H, mats.M, nmodes, 'smallestabs');
                lambda = diag(D);
                if any(~isfinite(lambda))
                    error('solve_bishop_ballooning_stability:NonFiniteFallback', ...
                        'smallestabs returned non-finite eigenvalues.');
                end
            catch ME2
                if mats.ni > 2500
                    error('solve_bishop_ballooning_stability:EigsFallbackFailed', ...
                        ['eigs fallback failed (%s), and the matrix is too large ', ...
                         'for a safe full eig fallback.'], ME2.message);
                end
                [Fi, D] = eig(full(mats.H), full(mats.M));
                lambda = diag(D);
            end
        end
    else
        [Fi, D] = eig(full(mats.H), full(mats.M));
        lambda = diag(D);
    end

    [~, ord] = sort(real(lambda), 'ascend');
    ord = ord(1:min(nmodes, numel(ord)));
    lambda = lambda(ord);
    Fi = Fi(:, ord);

    F = zeros(mats.n, numel(lambda));
    F(mats.interior, :) = Fi;

    for k = 1:size(F, 2)
        normM = sqrt(real(Fi(:,k)' * mats.M * Fi(:,k)));
        if normM > 0
            F(:,k) = F(:,k) ./ normM;
            Fi(:,k) = Fi(:,k) ./ normM;
        end

        [~, imax] = max(abs(F(:,k)));
        if F(imax,k) < 0
            F(:,k) = -F(:,k);
            Fi(:,k) = -Fi(:,k);
        end
    end

    sol = struct();
    sol.lambda = lambda;
    sol.Fi = Fi;
    sol.F = F;
end

%======================================================================
function target = choose_eigs_target(targetOpt, lambdaMinEstimate)
    targetString = lower(string(targetOpt));

    if targetString == "auto"
        delta = max(1e-4 * max(1, abs(lambdaMinEstimate)), 1e-8);
        target = lambdaMinEstimate - delta;
    else
        target = char(targetOpt);
    end
end

%======================================================================
function lambdaMin = estimate_lowest_eigenvalue(H, M)
    % H - sigma*M is positive definite exactly when sigma is below the
    % smallest generalized eigenvalue. Use sparse Cholesky tests to bracket
    % and bisect that threshold.
    if is_positive_definite(H)
        lo = 0;
        hi = 1;
        while is_positive_definite(H - hi.*M)
            hi = 2 * hi;
            if hi > 1e16
                error('solve_bishop_ballooning_stability:EigenBracketFailed', ...
                    'Could not bracket the lowest positive eigenvalue.');
            end
        end
    else
        hi = 0;
        lo = -1;
        while ~is_positive_definite(H - lo.*M)
            lo = 2 * lo;
            if abs(lo) > 1e16
                error('solve_bishop_ballooning_stability:EigenBracketFailed', ...
                    'Could not bracket the lowest negative eigenvalue.');
            end
        end
    end

    for it = 1:70
        mid = 0.5 * (lo + hi);
        if is_positive_definite(H - mid.*M)
            lo = mid;
        else
            hi = mid;
        end
    end

    lambdaMin = 0.5 * (lo + hi);
end

%======================================================================
function tf = is_positive_definite(A)
    A = 0.5 * (A + A.');
    [~, p] = chol(A);
    tf = (p == 0);
end

%======================================================================
function energies = compute_mode_energies(mats, sol)
    nm = numel(sol.lambda);

    bending = zeros(nm, 1);
    drive = zeros(nm, 1);
    deltaW = zeros(nm, 1);
    normM = zeros(nm, 1);

    for k = 1:nm
        fi = sol.Fi(:,k);
        bending(k) = real(fi' * mats.K * fi);
        drive(k) = real(fi' * mats.P * fi);
        deltaW(k) = real(fi' * mats.H * fi);
        normM(k) = real(fi' * mats.M * fi);
    end

    energies = struct();
    energies.bending = bending;
    energies.pressureCurvature = drive;
    energies.deltaW = deltaW;
    energies.norm = normM;
    energies.rayleigh = deltaW ./ normM;
end

%======================================================================
function plot_stability_solution(stab)
    figure('Color', 'w', 'Name', 'Bishop ballooning stability');
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(stab.theta, stab.F(:,1), 'LineWidth', 1.5);
    grid on;
    xlabel('theta');
    ylabel('F');
    title(sprintf('Lowest mode: lambda = %.6e', real(stab.lambda(1))));

    nexttile;
    plot(stab.theta, stab.mats.Acoef, 'LineWidth', 1.2);
    grid on;
    xlabel('theta');
    ylabel('A');
    title('Bending coefficient');

    nexttile;
    plot(stab.theta, stab.mats.Vcoef, 'LineWidth', 1.2);
    grid on;
    xlabel('theta');
    ylabel('V');
    title('Drive potential C/Bp');

    nexttile;
    bar(real(stab.lambda));
    grid on;
    xlabel('mode');
    ylabel('lambda');
    title(stab.status);
end
