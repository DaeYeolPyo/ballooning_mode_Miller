function bal_coef = evaluate_gcf(diag, theta_bal, opts)
    arguments
        diag (1,1) struct
        theta_bal (1,1) double
        opts.theta_bnd (1,1) double {mustBePositive} = 5*pi
        opts.B_N (1,1) double {mustBePositive} = 1.0
        opts.a_N (1,1) double {mustBePositive} = 1.0
    end

    required = {'theta','theta_psi','R','Z','Bp','B2','H', ...
        'normal_R','normal_Z','tangent_R','tangent_Z', ...
        'curvature_normal','dBp_dl','F','q','qprime','pprime'};
    missing = required(~isfield(diag, required));
    if ~isempty(missing)
        error('evaluate_gcf:MissingDiagnosticFields', ...
            'diag is missing required field(s): %s', strjoin(missing, ', '));
    end

    theta = diag.theta(:);
    if numel(theta) < 4 || any(diff(theta) <= 0)
        error('evaluate_gcf:BadThetaGrid', ...
            'diag.theta must be a strictly increasing closed periodic grid.');
    end

    period = theta(end) - theta(1);
    if abs(period - 2*pi) > 1.e-8
        error('evaluate_gcf:BadThetaPeriod', ...
            'diag.theta must span one complete 2*pi period.');
    end

    % Use a uniform extended ballooning grid with the same number of cells per
    % 2*pi period as the closed reference-surface grid.
    nperiod = numel(theta) - 1;
    dtheta = 2*pi/nperiod;
    nhalf = max(1, round(opts.theta_bnd/dtheta));
    theta_ext = (-nhalf:nhalf).'*dtheta;

    R = extend_periodic(theta, diag.R, theta_ext);
    Z = extend_periodic(theta, diag.Z, theta_ext);
    Bp = extend_periodic(theta, diag.Bp, theta_ext);
    B2_interpolated = extend_periodic(theta, diag.B2, theta_ext);
    H_interpolated = extend_periodic(theta, diag.H, theta_ext);
    theta_psi = extend_periodic(theta, diag.theta_psi, theta_ext);
    normal_R = extend_periodic(theta, diag.normal_R, theta_ext);
    normal_Z = extend_periodic(theta, diag.normal_Z, theta_ext);
    tangent_R = extend_periodic(theta, diag.tangent_R, theta_ext);
    tangent_Z = extend_periodic(theta, diag.tangent_Z, theta_ext);
    curvature_normal = extend_periodic( ...
        theta, diag.curvature_normal, theta_ext);
    dBp_dl = extend_periodic(theta, diag.dBp_dl, theta_ext);

    F = diag.F;
    q = diag.q;
    qprime = diag.qprime;
    pprime = diag.pprime;
    aN = opts.a_N;
    BN = opts.B_N;
    mu0 = 4*pi*1.e-7;

    Bphi = F./R;
    B2 = Bp.^2 + Bphi.^2;
    B = sqrt(B2);
    H = F./(q.*R.^2.*Bp);

    normal_norm = hypot(normal_R, normal_Z);
    tangent_norm = hypot(tangent_R, tangent_Z);
    normal_R = normal_R./normal_norm;
    normal_Z = normal_Z./normal_norm;
    tangent_R = tangent_R./tangent_norm;
    tangent_Z = tangent_Z./tangent_norm;
    delta_theta = theta_ext - theta_bal;
    I = qprime.*delta_theta + q.*theta_psi;

    % Clebsch-gradient norm in physical coordinates.
    K = 1./R.^2 + F^2./(R.^4.*Bp.^2) + (R.*Bp).^2.*I.^2;

    % P = mu0*p + B^2/2.  These normal and tangential derivatives follow
    % directly from the Mercier-Luc local equilibrium and force balance.
    P_normal = Bp.^2.*curvature_normal - F^2.*normal_R./R.^3;
    P_tangent = Bp.*dBp_dl - F^2.*tangent_R./R.^3;

    % This sign is tied to the repository convention B = grad(psi) x
    % grad(alpha), with alpha = zeta - q(theta-theta0).
    D = B2./(R.*Bp).*P_normal - F.*Bp.*I.*P_tangent;

    % Gaur normalization: grad_N = a_N*grad and
    % psi_N = psi/(a_N^2*B_N).
    dpbar_dpsiN = mu0*aN^2*pprime/BN;
    dq_dpsiN = aN^2*BN*qprime;
    BoverBN = B./BN;
    b_dot_gradN_theta = aN*F./(q.*R.^2.*B);
    gradN_alpha_sq = aN^2.*K;

    g = b_dot_gradN_theta.*gradN_alpha_sq./BoverBN;
    f = gradN_alpha_sq./(b_dot_gradN_theta.*BoverBN.^3);
    c = 2*aN*q*BN^2.*R.^2.*dpbar_dpsiN.*D ...
        ./ (F.*B2.^2);

    % Independent vector reconstruction in the local cylindrical orthonormal
    % basis [e_R,e_phi,e_Z].  These checks do not reuse K or D.
    normal = [normal_R, zeros(size(R)), normal_Z];
    tangent = [tangent_R, zeros(size(R)), tangent_Z];
    grad_psi = (R.*Bp).*normal;
    grad_theta = (R.*Bp.*theta_psi).*normal + H.*tangent;
    grad_zeta = [zeros(size(R)), 1./R, zeros(size(R))];
    grad_alpha = grad_zeta - q.*grad_theta ...
        - (qprime.*delta_theta).*grad_psi;

    Bvec = [Bp.*tangent_R, Bphi, Bp.*tangent_Z];
    bvec = Bvec./B;
    grad_P = P_normal.*normal + P_tangent.*tangent;
    D_vector = sum(cross(Bvec, grad_P, 2).*grad_alpha, 2);

    gradN_alpha = aN.*grad_alpha;
    gradN_H = (aN/BN^2).*grad_P;
    gradN_alpha_sq_vector = sum(gradN_alpha.^2, 2);
    b_dot_gradN_theta_vector = sum(bvec.*(aN.*grad_theta), 2);
    cross_term = sum(cross(bvec, gradN_H, 2).*gradN_alpha, 2);

    g_vector = b_dot_gradN_theta_vector.*gradN_alpha_sq_vector./BoverBN;
    f_vector = gradN_alpha_sq_vector ...
        ./(b_dot_gradN_theta_vector.*BoverBN.^3);
    c_vector = 2./b_dot_gradN_theta_vector./BoverBN.^4 ...
        .*dpbar_dpsiN.*cross_term;

    validation = struct();
    validation.basis_normal_norm_error = max(abs( ...
        sum(normal.^2, 2) - 1));
    validation.basis_tangent_norm_error = max(abs( ...
        sum(tangent.^2, 2) - 1));
    validation.basis_orthogonality_error = max(abs( ...
        sum(normal.*tangent, 2)));
    validation.B2_interpolation_relative_error = relative_error( ...
        B2_interpolated, B2);
    validation.H_interpolation_relative_error = relative_error( ...
        H_interpolated, H);
    validation.K_relative_error = relative_error(K, sum(grad_alpha.^2, 2));
    validation.D_relative_error = relative_error(D, D_vector);
    validation.b_dot_gradN_theta_relative_error = relative_error( ...
        b_dot_gradN_theta, b_dot_gradN_theta_vector);
    validation.g_relative_error = relative_error(g, g_vector);
    validation.f_relative_error = relative_error(f, f_vector);
    validation.c_relative_error = relative_error(c, c_vector);
    validation.all_finite = all(isfinite([g; c; f; I; K; D]));

    bal_coef = struct();
    bal_coef.model = 'miller-mercier-luc';
    bal_coef.theta = theta_ext;
    bal_coef.theta0 = theta_bal;
    bal_coef.theta_bnd_requested = opts.theta_bnd;
    bal_coef.theta_bnd_actual = max(abs(theta_ext));
    bal_coef.g = g;
    bal_coef.c = c;
    bal_coef.f = f;
    bal_coef.I = I;
    bal_coef.K = K;
    bal_coef.D = D;
    bal_coef.P_normal = P_normal;
    bal_coef.P_tangent = P_tangent;

    bal_coef.R = R;
    bal_coef.Z = Z;
    bal_coef.Bp = Bp;
    bal_coef.Bphi = Bphi;
    bal_coef.B2 = B2;
    bal_coef.B = B;
    bal_coef.BoverBN = BoverBN;
    bal_coef.H = H;
    bal_coef.theta_psi = theta_psi;
    bal_coef.normal = normal;
    bal_coef.tangent = tangent;
    bal_coef.curvature_normal = curvature_normal;
    bal_coef.dBp_dl = dBp_dl;

    bal_coef.F = F;
    bal_coef.q = q;
    bal_coef.qprime = qprime;
    bal_coef.pprime = pprime;
    bal_coef.FFprime = diag.FFprime;
    bal_coef.aN = aN;
    bal_coef.BN = BN;
    bal_coef.dq_dpsiN = dq_dpsiN;
    bal_coef.dpbar_dpsiN = dpbar_dpsiN;

    bal_coef.grad_psi = grad_psi;
    bal_coef.grad_theta = grad_theta;
    bal_coef.grad_zeta = grad_zeta;
    bal_coef.grad_alpha = grad_alpha;
    bal_coef.gradN_alpha = gradN_alpha;
    bal_coef.gradN_alpha_sq = gradN_alpha_sq;
    bal_coef.b_dot_gradN_theta = b_dot_gradN_theta;
    bal_coef.Bvec = Bvec;
    bal_coef.bvec = bvec;
    bal_coef.grad_P = grad_P;
    bal_coef.gradN_H = gradN_H;
    bal_coef.cross_term = cross_term;
    bal_coef.D_vector = D_vector;
    bal_coef.g_vector = g_vector;
    bal_coef.c_vector = c_vector;
    bal_coef.f_vector = f_vector;
    bal_coef.validation = validation;
end

function Q_ext = extend_periodic(theta, Q, theta_ext)
    theta = theta(:);
    if size(Q, 1) ~= numel(theta)
        Q = Q.';
    end
    if size(Q, 1) ~= numel(theta)
        error('evaluate_gcf:BadPeriodicField', ...
            'A periodic field must have one row per theta sample.');
    end

    period = theta(end) - theta(1);
    theta_core = theta(1:end-1);
    Q_core = Q(1:end-1, :);
    theta_interp = [theta_core - period; theta_core; theta_core + period];
    Q_interp = [Q_core; Q_core; Q_core];
    theta_query = mod(theta_ext - theta(1), period) + theta(1);
    Q_ext = interp1(theta_interp, Q_interp, theta_query, 'pchip');
end

function err = relative_error(actual, expected)
    scale = max(max(abs(expected(:))), eps);
    err = max(abs(actual(:) - expected(:)))/scale;
end
