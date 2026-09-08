function coeff = calculate_ballooning_coeffs(map, metrics, target_psiN, opts)
%CALCULATE_BALLOONING_COEFFS Compute PEST g, c, and f on one flux surface.
%
% The returned coefficients satisfy
%
%   d/dtheta (g dX/dtheta) + c X = lambda f X.
%
% Only the requested flux surface is returned.  Neighboring surfaces in
% map and metrics are nevertheless required to evaluate dq/dpsi and
% d(B^2)/dpsi.  The target surface must therefore be present in map.psiN
% and must have at least two neighboring surfaces on either side.
%
% Inputs
%   map          : maps.PEST from build_SFL_maps
%   metrics      : compute_metrics(map, eqfunc) output
%   target_psiN  : normalized poloidal-flux label of the target surface
%
% Name-value options
%   Theta0       : ballooning angle theta_0 (default 0)
%   ThetaExt     : strictly increasing evaluation grid.  If empty, one
%                  complete period [theta(1), theta(1)+2*pi] is returned.
%
% Output
%   coeff.theta, coeff.g, coeff.c, coeff.f are column vectors ready for
%   construct_matrix.  coeff.base contains the periodic building blocks.

    arguments
        map (1,1) struct
        metrics (1,1) struct
        target_psiN (1,1) double
        opts.Theta0 (1,1) double = 0
        opts.ThetaExt double = []
        opts.use_salpha (1,1) logical = false
        opts.s (1,1) double = 0
        opts.alpha (1, 1) double = 0
    end

    validate_inputs(map, metrics);

    psiN = map.psiN(:);
    psi = map.psi(:);
    thetaBase = map.theta(:).';

    [~, k] = min(abs(psiN - target_psiN));
    targetTolerance = 100 * eps(max(1, max(abs(psiN))));

    if abs(psiN(k) - target_psiN) > targetTolerance
        error('calculate_ballooning_coeffs:TargetNotOnGrid', ...
            ['target_psiN = %.16g is not present in map.psiN. ', ...
             'Include the target surface when constructing the SFL map. ', ...
             'The nearest available surface is %.16g.'], ...
            target_psiN, psiN(k));
    end

    Nrho = numel(psiN);
    if Nrho < 5 || k < 3 || k > Nrho - 2
        error('calculate_ballooning_coeffs:InsufficientRadialStencil', ...
            ['The target surface needs two neighboring flux surfaces on ', ...
             'each side for the five-point local cubic radial fit.']);
    end

    if any(diff(psiN) <= 0)
        error('calculate_ballooning_coeffs:BadPsiNGrid', ...
            'map.psiN must be strictly increasing.');
    end

    if any(diff(psi) == 0) || ...
            ~(all(diff(psi) > 0) || all(diff(psi) < 0))
        error('calculate_ballooning_coeffs:BadPsiGrid', ...
            'map.psi must be strictly monotonic.');
    end

    if isempty(opts.ThetaExt)
        thetaExt = [thetaBase, thetaBase(1) + 2*pi];
    else
        if ~isvector(opts.ThetaExt)
            error('calculate_ballooning_coeffs:BadThetaExt', ...
                'ThetaExt must be a vector.');
        end
        thetaExt = opts.ThetaExt(:).';
    end

    if numel(thetaExt) < 2 || any(~isfinite(thetaExt)) || ...
            any(diff(thetaExt) <= 0)
        error('calculate_ballooning_coeffs:BadThetaExt', ...
            'ThetaExt must be finite and strictly increasing.');
    end

    % The q used in alpha must be a smooth equilibrium profile.  Use a
    % separately stored profile when available; otherwise fall back to
    % map.q (currently ang.q_used in build_SFL_maps).
    if isfield(map, 'q_profile')
        qProfile = map.q_profile(:);
        qSource = "profile";
    else
        qProfile = map.q(:);
        qSource = "map.q";
    end

    pprimeProfile = map.pprime(:);

    if numel(qProfile) ~= Nrho || numel(pprimeProfile) ~= Nrho
        error('calculate_ballooning_coeffs:BadProfileSize', ...
            'q and pprime must contain one value per flux surface.');
    end

    % q' and B2_psi are derivatives with respect to map.psi, not psiN.
    q = qProfile(k);

    if opts.use_salpha
        
    else
        qprime = local_cubic_derivative(psi, qProfile, k);
        pprime = pprimeProfile(k);
    end

    B2all = metrics.B2;
    B2 = B2all(k, :);
    B2_psi = local_cubic_derivative(psi, B2all, k);
    B2_theta = periodic_theta_derivative(B2, thetaBase);

    J = metrics.Jacobian(k, :);
    R = map.R(k, :);

    if any(~isfinite(J)) || any(J <= 0)
        error('calculate_ballooning_coeffs:NonPositiveJacobian', ...
            ['The PEST Jacobian must be finite and positive. ', ...
             'Fix the flux-surface orientation instead of silently ', ...
             'replacing it by abs(J).']);
    end

    if any(~isfinite(B2)) || any(B2 <= 0)
        error('calculate_ballooning_coeffs:BadB2', ...
            'B2 must be finite and positive on the target surface.');
    end

    gcon = metrics.g_contra;
    gcov = metrics.g_cov;

    % PEST: alpha_i = (-q' DeltaTheta, -q, 1).  Because
    % g^{psi,zeta}=g^{theta,zeta}=0, K=|grad alpha|^2 has the form
    % K0 + DeltaTheta*K1 + DeltaTheta^2*K2.
    K0 = q^2 .* gcon.thetatheta(k, :) + gcon.zetazeta(k, :);
    K1 = 2*q*qprime .* gcon.psitheta(k, :);
    K2 = qprime^2 .* gcon.psipsi(k, :);

    denomGF = J .* B2;

    g0 = K0 ./ denomGF;
    g1 = K1 ./ denomGF;
    g2 = K2 ./ denomGF;

    f0 = J .* K0 ./ B2;
    f1 = J .* K1 ./ B2;
    f2 = J .* K2 ./ B2;

    % Pressure-curvature coefficient in axisymmetric PEST coordinates.
    Gpar = gcov.thetatheta(k, :) + q^2 .* R.^2;
    denomC = J .* B2.^2;

    c0 = pprime ./ denomC .* ( ...
          Gpar .* (B2_psi + 2*pprime) ...
        - gcov.psitheta(k, :) .* B2_theta);

    c1 = -pprime ./ denomC .* ...
        (q*qprime .* R.^2 .* B2_theta);

    % Periodic geometry is evaluated on thetaExt modulo 2*pi.  Only the
    % explicit DeltaTheta terms remain non-periodic.
    deltaTheta = thetaExt - opts.Theta0;

    g0e = periodic_eval(thetaBase, g0, thetaExt);
    g1e = periodic_eval(thetaBase, g1, thetaExt);
    g2e = periodic_eval(thetaBase, g2, thetaExt);

    f0e = periodic_eval(thetaBase, f0, thetaExt);
    f1e = periodic_eval(thetaBase, f1, thetaExt);
    f2e = periodic_eval(thetaBase, f2, thetaExt);

    c0e = periodic_eval(thetaBase, c0, thetaExt);
    c1e = periodic_eval(thetaBase, c1, thetaExt);

    K0e = periodic_eval(thetaBase, K0, thetaExt);
    K1e = periodic_eval(thetaBase, K1, thetaExt);
    K2e = periodic_eval(thetaBase, K2, thetaExt);

    K = K0e + deltaTheta .* K1e + deltaTheta.^2 .* K2e;
    g = g0e + deltaTheta .* g1e + deltaTheta.^2 .* g2e;
    f = f0e + deltaTheta .* f1e + deltaTheta.^2 .* f2e;
    c = c0e + deltaTheta .* c1e;

    if any(~isfinite(K)) || any(K <= 0) || ...
            any(~isfinite(g)) || any(g <= 0) || ...
            any(~isfinite(f)) || any(f <= 0) || ...
            any(~isfinite(c))
        error('calculate_ballooning_coeffs:InvalidCoefficient', ...
            ['Non-finite or non-positive K, g, or f was produced. ', ...
             'Check the metric, Jacobian, q'', and radial derivatives.']);
    end

    coeff = struct();

    coeff.target_psiN = psiN(k);
    coeff.target_psi = psi(k);
    coeff.target_index = k;
    coeff.theta0 = opts.Theta0;

    % Column-vector outputs can be passed directly to construct_matrix.
    coeff.theta = thetaExt(:);
    coeff.g = g(:);
    coeff.c = c(:);
    coeff.f = f(:);
    coeff.K = K(:);

    coeff.q = q;
    coeff.qprime = qprime;
    coeff.q_source = qSource;
    coeff.pprime = pprime;

    coeff.base = struct();
    coeff.base.theta = thetaBase(:);
    coeff.base.B2 = B2(:);
    coeff.base.B2_psi = B2_psi(:);
    coeff.base.B2_theta = B2_theta(:);

    coeff.base.K0 = K0(:);
    coeff.base.K1 = K1(:);
    coeff.base.K2 = K2(:);

    coeff.base.g0 = g0(:);
    coeff.base.g1 = g1(:);
    coeff.base.g2 = g2(:);

    coeff.base.f0 = f0(:);
    coeff.base.f1 = f1(:);
    coeff.base.f2 = f2(:);

    coeff.base.c0 = c0(:);
    coeff.base.c1 = c1(:);

    coeff.diagnostics = struct();
    coeff.diagnostics.minK = min(K);
    coeff.diagnostics.minG = min(g);
    coeff.diagnostics.minF = min(f);
    coeff.diagnostics.maxAbsC = max(abs(c));

    if isfield(map, 'lambda')
        coeff.diagnostics.maxAbsLambda = max(abs(map.lambda(k, :)));
    end

    if isfield(metrics, 'B2_metric')
        B2metric = metrics.B2_metric(k, :);
        coeff.diagnostics.B2MetricRelativeRms = ...
            rms(B2metric - B2) / max(rms(B2), eps);
    end
end

function validate_inputs(map, metrics)
    mapFields = {'psi', 'psiN', 'theta', 'R', 'q', 'pprime'};
    metricFields = {'B2', 'Jacobian', 'g_cov', 'g_contra'};

    for k = 1:numel(mapFields)
        if ~isfield(map, mapFields{k})
            error('calculate_ballooning_coeffs:MissingMapField', ...
                'map.%s is required.', mapFields{k});
        end
    end

    for k = 1:numel(metricFields)
        if ~isfield(metrics, metricFields{k})
            error('calculate_ballooning_coeffs:MissingMetricField', ...
                'metrics.%s is required.', metricFields{k});
        end
    end

    [Nrho, Ntheta] = size(map.R);

    if ~isequal(size(metrics.B2), [Nrho, Ntheta]) || ...
            ~isequal(size(metrics.Jacobian), [Nrho, Ntheta])
        error('calculate_ballooning_coeffs:GridSizeMismatch', ...
            'map and metrics must use the same Nrho-by-Ntheta grid.');
    end

    if numel(map.psi) ~= Nrho || numel(map.psiN) ~= Nrho || ...
            numel(map.theta) ~= Ntheta
        error('calculate_ballooning_coeffs:GridSizeMismatch', ...
            'psi, psiN, and theta sizes do not match map.R.');
    end

    requiredCov = {'psitheta', 'thetatheta'};
    requiredCon = {'psipsi', 'psitheta', 'thetatheta', 'zetazeta'};

    for k = 1:numel(requiredCov)
        name = requiredCov{k};
        if ~isfield(metrics.g_cov, name) || ...
                ~isequal(size(metrics.g_cov.(name)), [Nrho, Ntheta])
            error('calculate_ballooning_coeffs:MissingMetricComponent', ...
                'metrics.g_cov.%s is missing or has the wrong size.', name);
        end
    end

    for k = 1:numel(requiredCon)
        name = requiredCon{k};
        if ~isfield(metrics.g_contra, name) || ...
                ~isequal(size(metrics.g_contra.(name)), [Nrho, Ntheta])
            error('calculate_ballooning_coeffs:MissingMetricComponent', ...
                'metrics.g_contra.%s is missing or has the wrong size.', name);
        end
    end
end

function dydx = local_cubic_derivative(x, y, centerIndex)
%LOCAL_CUBIC_DERIVATIVE Five-point least-squares cubic derivative.
    x = x(:);

    if size(y, 1) ~= numel(x)
        error('calculate_ballooning_coeffs:BadDerivativeInput', ...
            'The first dimension of y must match x.');
    end

    idx = centerIndex + (-2:2);
    dx = x(idx) - x(centerIndex);

    % The coefficient of dx is the derivative at the stencil center.
    V = [ones(5,1), dx, dx.^2, dx.^3];
    coefficient = V \ y(idx, :);
    dydx = coefficient(2, :);
end

function dA = periodic_theta_derivative(A, theta)
%PERIODIC_THETA_DERIVATIVE Fourier derivative on a uniform 2*pi grid.
    A = A(:).';
    theta = theta(:).';

    N = numel(theta);
    if numel(A) ~= N || N < 4
        error('calculate_ballooning_coeffs:BadThetaDerivativeInput', ...
            'A and theta must have the same length of at least four.');
    end

    dtheta = diff([theta, theta(1) + 2*pi]);
    if max(abs(dtheta - mean(dtheta))) > ...
            1e-10 * max(1, abs(mean(dtheta)))
        error('calculate_ballooning_coeffs:NonUniformTheta', ...
            'The PEST theta grid must be uniform for Fourier differentiation.');
    end

    if mod(N, 2) == 0
        modeNumber = [0:N/2-1, 0, -N/2+1:-1];
    else
        modeNumber = [0:(N-1)/2, -(N-1)/2:-1];
    end

    dA = real(ifft(1i .* modeNumber .* fft(A)));
end

function yq = periodic_eval(theta, y, thetaq)
%PERIODIC_EVAL Periodic PCHIP evaluation of a one-period scalar field.
    theta = theta(:).';
    y = y(:).';
    thetaq = thetaq(:).';

    theta0 = theta(1);
    thetaWrapped = mod(thetaq - theta0, 2*pi) + theta0;

    thetaExtended = [theta, theta0 + 2*pi];
    yExtended = [y, y(1)];

    yq = interp1(thetaExtended, yExtended, thetaWrapped, 'pchip');
end
