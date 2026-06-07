function sol = solve_ballooning_eigenvalue(bal, varargin)
%SOLVE_BALLOONING_EIGENVALUE Solve ideal ballooning eigenvalue problem.
%
%   sol = solve_ballooning_eigenvalue(bal)
%   sol = solve_ballooning_eigenvalue(bal,'Name',value,...)
%
% Input
%   bal : output of eq_ballooning_coefficients.m. Newer bal structs contain
%         one-period geometry fields, from which g,c,f are rebuilt on the
%         extended ballooning coordinate. Older structs with only bal.g,c,f
%         are still supported by periodic interpolation.
%
% Solves
%   d/dtheta ( g dX/dtheta ) + c X = lambda f X
%
% on theta_b in [-ThetaB, ThetaB] with X(-ThetaB)=X(ThetaB)=0.
%
% Periodic equilibrium geometry is extended along the field line, while the
% magnetic-shear part of grad(alpha_t) is kept non-periodic.
%
% Name-value options
%   ThetaB        : ballooning domain half-width, default 5*pi
%   N             : number of total grid points including boundaries, default 501
%                   Use odd N as in the paper. Interior unknowns are N-2.
%   Theta0        : ballooning parameter theta0 in alpha_t, default bal.theta0
%   AlphaShift    : optional shift of the equilibrium theta origin, default 0
%   NEigs         : number of eigenvalues requested, default 6
%   Which         : eigs selector, default 'largestreal'
%   UseSparse     : use sparse matrix, default true
%   FullEigThreshold : use full eig when interior size is below this, default 800
%   RefineRayleigh: recompute lambda with the assembled discrete Rayleigh
%                   quotient, default true
%   Plot          : true/false, default false
%
% Output
%   sol.lambda          : selected eigenvalue from matrix solve
%   sol.lambda_refined  : Rayleigh quotient value using returned eigenfunction
%   sol.theta_b         : full ballooning grid including boundaries
%   sol.X               : eigenfunction including zero boundaries
%   sol.A               : inv(M)*K operator consistent with (3.3)/(3.8)
%   sol.A_appendix      : -sol.A, the sign convention printed in Appendix A
%   sol.theta_int       : interior theta grid
%   sol.g, sol.c, sol.f : coefficients on full grid
%   sol.g_half          : g at half points
%
% Notes
%   Matrix convention used here:
%       K X = lambda M X
%   where K discretizes d/dtheta(g dX/dtheta) + c X and M=diag(f).
%   This convention is consistent with the Rayleigh quotient
%       lambda = int(c|X|^2 - g|X'|^2)/int(f|X|^2).

    defaultTheta0 = 0;
    if isfield(bal, 'theta0') && isscalar(bal.theta0) && isfinite(bal.theta0)
        defaultTheta0 = bal.theta0;
    end

    p = inputParser;
    addParameter(p, 'ThetaB', 5*pi, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'N', 501, @(x)isnumeric(x)&&isscalar(x)&&x>=5);
    addParameter(p, 'Theta0', defaultTheta0, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'AlphaShift', 0, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'NEigs', 6, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
    addParameter(p, 'Which', 'largestreal', @(x)ischar(x)||isstring(x));
    addParameter(p, 'UseSparse', true, @(x)islogical(x)&&isscalar(x));
    addParameter(p, 'FullEigThreshold', 800, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
    addParameter(p, 'RefineRayleigh', true, @(x)islogical(x)&&isscalar(x));
    addParameter(p, 'Plot', false, @(x)islogical(x)&&isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    N = round(opt.N);
    if mod(N,2) == 0
        warning('solve_ballooning_eigenvalue:EvenN', ...
            'N should be odd. Changing N=%d to N=%d.', N, N+1);
        N = N + 1;
    end

    theta_b = linspace(-opt.ThetaB, opt.ThetaB, N).';
    dth = theta_b(2) - theta_b(1);

    theta_eval = theta_b + opt.AlphaShift;
    theta_half = 0.5*(theta_b(1:end-1) + theta_b(2:end)) + opt.AlphaShift;

    % Rebuild the coefficients on the extended ballooning coordinate.
    % The equilibrium geometry is periodic in theta, but grad(alpha_t)
    % contains the secular magnetic-shear term -q'(theta-theta0) grad(psi).
    % Therefore g,c,f must not be obtained by periodically repeating a
    % single-period g,c,f table.
    if has_ballooning_geometry(bal)
        [g, c, f] = evaluate_ballooning_coefficients(bal, theta_eval, opt.Theta0);
        [g_half, ~, ~] = evaluate_ballooning_coefficients(bal, theta_half, opt.Theta0);
    else
        warning('solve_ballooning_eigenvalue:PeriodicCoefficientFallback', ...
            ['bal lacks geometry fields. Falling back to periodic g,c,f ', ...
             'extension, which omits the q''*(theta-theta0) shear growth.']);
        [g, c, f] = interpolate_periodic_coefficients(bal.theta(:), bal.g(:), bal.c(:), bal.f(:), ...
                                                       theta_eval);
        [g_half, ~, ~] = interpolate_periodic_coefficients(bal.theta(:), bal.g(:), bal.c(:), bal.f(:), ...
                                                           theta_half);
    end

    % Unknowns are interior full-grid points: j = 2,...,N-1 in MATLAB indexing.
    M = N - 2;
    theta_int = theta_b(2:end-1);
    c_int = c(2:end-1);
    f_int = f(2:end-1);

    lowerK = zeros(M,1);
    diagK  = zeros(M,1);
    upperK = zeros(M,1);

    for m = 1:M
        % Physical full-grid index j = m+1.
        % g_minus = g_{j-1/2}; g_plus = g_{j+1/2}.
        g_minus = g_half(m);
        g_plus  = g_half(m+1);
        cj = c_int(m);

        % Discretization of
        %   d/dtheta(g dX/dtheta) + c X = lambda f X
        % gives
        %   K(j,j-1)= g_minus/dtheta^2
        %   K(j,j)  = c_j - (g_minus+g_plus)/dtheta^2
        %   K(j,j+1)= g_plus/dtheta^2
        diagK(m) = cj - (g_minus + g_plus)/dth^2;
        if m > 1
            lowerK(m) = g_minus/dth^2;
        end
        if m < M
            upperK(m) = g_plus/dth^2;
        end
    end

    if opt.UseSparse
        K = spdiags([[lowerK(2:end);0], diagK, [0;upperK(1:end-1)]], -1:1, M, M);
        Mmat = spdiags(f_int, 0, M, M);
    else
        K = diag(diagK) + diag(upperK(1:end-1),1) + diag(lowerK(2:end),-1);
        Mmat = diag(f_int);
    end

    % For backward compatibility, A stores the divided operator inv(M)*K.
    % Appendix A prints the same tridiagonal stencil with the opposite sign.
    A = Mmat \ K;

    % eigs option string differs across MATLAB versions. Try robustly.
    neigs = min(round(opt.NEigs), M-2);
    which = char(opt.Which);

    usedFullEig = M <= opt.FullEigThreshold;
    if usedFullEig
        [V,D] = eig(full(K), full(Mmat));
        evals = diag(D);
    else
        try
            [V,D] = eigs(K, Mmat, neigs, which);
            evals = diag(D);
            if ~any(isfinite(evals))
                error('solve_ballooning_eigenvalue:EigsNoFiniteValues', ...
                    'eigs returned no finite eigenvalues.');
            end
        catch
            % Fallback: full eigenvalue solve for moderate matrix sizes.
            [V,D] = eig(full(K), full(Mmat));
            evals = diag(D);
            usedFullEig = true;
        end
    end

    % Pick largest real part eigenvalue.
    finite = isfinite(evals);
    if ~any(finite)
        error('solve_ballooning_eigenvalue:NoFiniteEigenvalues', ...
            'No finite eigenvalues were found.');
    end
    evals_select = evals;
    evals_select(~finite) = -Inf;
    [~, imax] = max(real(evals_select));
    lambda = evals(imax);
    Xint = V(:,imax);

    % Normalize sign/scale.
    X = [0; Xint; 0];
    [~, im] = max(abs(X));
    X = X / X(im);
    if real(X(im)) < 0
        X = -X;
    end
    X = real(X);

    lambda_refined = NaN;
    if opt.RefineRayleigh
        lambda_refined = discrete_rayleigh_refine(X(2:end-1), K, Mmat);
    end

    sol = struct();
    sol.lambda = lambda;
    sol.lambda_refined = lambda_refined;
    sol.theta_b = theta_b;
    sol.theta_int = theta_int;
    sol.X = X;
    sol.A = A;
    sol.A_appendix = -A;
    sol.K = K;
    sol.M = Mmat;
    sol.g = g;
    sol.c = c;
    sol.f = f;
    sol.g_half = g_half;
    sol.dtheta = dth;
    sol.theta0 = opt.Theta0;
    sol.theta_eval = theta_eval;
    sol.used_full_eig = usedFullEig;
    sol.options = opt;

    if opt.Plot
        figure;
        plot(theta_b, X, 'LineWidth', 1.5);
        grid on;
        xlabel('\theta_b');
        ylabel('X');
        title(sprintf('Ballooning eigenfunction, lambda = %.6g, refined = %.6g', ...
            real(lambda), lambda_refined));
    end
end

%==========================================================================
function tf = has_ballooning_geometry(bal)
    req = {'theta','q','dq_dpsin','dpbar_dpsin','BoverBN', ...
           'b_dot_gradN_theta','bvec','gradN_H', ...
           'gradN_psin','gradN_theta','gradN_zeta'};
    tf = all(isfield(bal, req));
end

%==========================================================================
function [g, c, f] = evaluate_ballooning_coefficients(bal, theta_eval, theta0)
    theta_eval = theta_eval(:);
    theta_base = bal.theta(:);

    BoverBN = interp_periodic(theta_base, bal.BoverBN(:), theta_eval);
    b_dot_gradN_theta = interp_periodic(theta_base, bal.b_dot_gradN_theta(:), theta_eval);
    bvec = interp_periodic(theta_base, bal.bvec, theta_eval);
    gradN_H = interp_periodic(theta_base, bal.gradN_H, theta_eval);
    gradN_psin = interp_periodic(theta_base, bal.gradN_psin, theta_eval);
    gradN_theta = interp_periodic(theta_base, bal.gradN_theta, theta_eval);
    gradN_zeta = interp_periodic(theta_base, bal.gradN_zeta, theta_eval);

    shear_term = bal.dq_dpsin .* (theta_eval - theta0);
    if isfield(bal, 'local_shear_model') && strcmpi(string(bal.local_shear_model), "miller-s-alpha")
        % Miller/Greene-Chance local equilibria contain a pressure-gradient
        % contribution to the integrated local shear. In the circular
        % large-aspect-ratio limit this gives s*theta - alpha*sin(theta),
        % not just s*theta. The Miller adapter stores q*alpha as the
        % periodic pressure-shear amplitude.
        shear_term = shear_term ...
            - interp_periodic(theta_base, bal.local_shear_periodic_alpha(:) .* sin(theta_base), theta_eval);
    end

    gradN_alpha = gradN_zeta ...
        - bal.q .* gradN_theta ...
        - shear_term .* gradN_psin;

    gradN_alpha_sq = sum(gradN_alpha.^2, 2);
    cross_term = sum(cross(bvec, gradN_H, 2) .* gradN_alpha, 2);

    g = b_dot_gradN_theta .* gradN_alpha_sq ./ BoverBN;
    f = gradN_alpha_sq ./ (b_dot_gradN_theta .* BoverBN.^3);
    c = 2 ./ b_dot_gradN_theta ./ BoverBN.^4 ...
        .* bal.dpbar_dpsin .* cross_term;
end

%==========================================================================
function [gq, cq, fq] = interpolate_periodic_coefficients(theta, g, c, f, thetaq)
    gq = interp_periodic(theta, g(:), thetaq);
    cq = interp_periodic(theta, c(:), thetaq);
    fq = interp_periodic(theta, f(:), thetaq);
end

%==========================================================================
function vq = interp_periodic(theta, v, thetaq)
    theta = theta(:);
    if size(v, 1) ~= numel(theta)
        v = v.';
    end
    if size(v, 1) ~= numel(theta)
        error('interp_periodic:BadSize', ...
            'Value table must have one row per theta sample.');
    end

    th = mod(theta, 2*pi);
    [th, idx] = sort(th);
    v = v(idx, :);
    [th, v] = unique_periodic_samples(th, v);

    th_ext = [th-2*pi; th; th+2*pi];
    v_ext = [v; v; v];
    tq = mod(thetaq(:), 2*pi);

    vq = interp1(th_ext, v_ext, tq, 'pchip');
end

%==========================================================================
function [th_u, v_u] = unique_periodic_samples(th, v)
    tol = 1e-12;
    if isempty(th)
        th_u = th;
        v_u = v;
        return;
    end
    group = cumsum([1; abs(diff(th)) > tol]);
    n = group(end);
    th_u = zeros(n,1);
    v_u = zeros(n, size(v,2));
    for k = 1:n
        m = group == k;
        th_u(k) = mean(th(m));
        v_u(k,:) = mean(v(m,:), 1, 'omitnan');
    end

    % Remove possible 2pi duplicate of 0.
    if numel(th_u) > 1 && abs((th_u(end)-th_u(1)) - 2*pi) < 1e-10
        th_u(end) = [];
        v_u(end,:) = [];
    end
end

%==========================================================================
function lam = discrete_rayleigh_refine(Xint, K, Mmat)
    num = Xint' * (K * Xint);
    den = Xint' * (Mmat * Xint);
    lam = real(num / den);
end
