function G = compute_metrics(map, eqfunc)
%COMPUTE_METRICS Metric, B components, and J components on SFL grid.
    mu0 = 4*pi*1e-7;

    psi   = map.psi(:);
    theta = map.theta(:).';

    R = map.R;
    Z = map.Z;
    lam = map.lambda;

    [Nrho, Ntheta] = size(R);

    if ~isfield(map, 'F')
        error('compute_metrics:MissingF', ...
              'map.F is required for direct B2 calculation.');
    end
    
    if numel(map.F) ~= Nrho
        error('compute_metrics:BadFSize', ...
              'map.F must have one value per flux surface.');
    end
    
    psiR = eqfunc.psiR(R, Z);
    psiZ = eqfunc.psiZ(R, Z);
    
    gradPsi2 = psiR.^2 + psiZ.^2;
    
    F2d = repmat(map.F(:), 1, Ntheta);

    q = map.q(:);
    q2d = repmat(q, 1, numel(theta));

    R_psi = d_dpsi(R, psi);
    Z_psi = d_dpsi(Z, psi);
    l_psi = d_dpsi(lam, psi);

    R_th = d_dtheta(R, theta);
    Z_th = d_dtheta(Z, theta);
    l_th = d_dtheta(lam, theta);

    % Positive plasma-physics Jacobian:
    % dV = sqrtg dpsi dtheta dzeta
    Jacobian = R .* (R_psi .* Z_th - Z_psi .* R_th);

    if median(Jacobian(:), 'omitnan') < 0
        warning('compute_metric_B_J:NegativeJacobian', ...
                'Jacobian is mostly negative. Check theta orientation.');
    end

    Jacobian_abs = abs(Jacobian);

    % Covariant metric g_ij = e_i dot e_j.
    g11 = R_psi.^2 + Z_psi.^2 + R.^2 .* l_psi.^2;
    g12 = R_psi.*R_th + Z_psi.*Z_th + R.^2 .* l_psi.*l_th;
    g13 = -R.^2 .* l_psi;

    g22 = R_th.^2 + Z_th.^2 + R.^2 .* l_th.^2;
    g23 = -R.^2 .* l_th;

    g33 = R.^2;

    % Contravariant metric g^ij = inverse(g_ij).
    detg = g11 .* (g22 .* g33 - g23.^2) ...
         - g12 .* (g12 .* g33 - g13 .* g23) ...
         + g13 .* (g12 .* g23 - g13 .* g22);

    if any(~isfinite(detg(:))) || any(abs(detg(:)) < eps)
        error('compute_metric_B_J:SingularMetric', ...
              'Metric tensor is singular or invalid at some grid points.');
    end

    gcon11 = (g22 .* g33 - g23.^2) ./ detg;
    gcon12 = (g13 .* g23 - g12 .* g33) ./ detg;
    gcon13 = (g12 .* g23 - g13 .* g22) ./ detg;

    gcon22 = (g11 .* g33 - g13.^2) ./ detg;
    gcon23 = (g12 .* g13 - g11 .* g23) ./ detg;

    gcon33 = (g11 .* g22 - g12.^2) ./ detg;

    % Magnetic field contravariant components in straight-field-line coords.
    Bcon.psi   = zeros(size(R));
    Bcon.theta = 1 ./ Jacobian_abs;
    Bcon.zeta  = q2d ./ Jacobian_abs;

    % Magnetic field covariant components B_i = g_ij B^j.
    Bcov.psi   = g12 .* Bcon.theta + g13 .* Bcon.zeta;
    Bcov.theta = g22 .* Bcon.theta + g23 .* Bcon.zeta;
    Bcov.zeta  = g23 .* Bcon.theta + g33 .* Bcon.zeta;

    Bp2   = gradPsi2 ./ R.^2;
    Bphi2 = F2d.^2   ./ R.^2;

    B2_direct = Bp2 + Bphi2;

    B2_metric = Bcov.psi .* Bcon.psi + ...
         Bcov.theta .* Bcon.theta + ...
         Bcov.zeta .* Bcon.zeta;

    % Current density from curl B:
    % mu0 J^i = (curl B)^i.
    dBpsi_dth = d_dtheta(Bcov.psi, theta);
    dBtheta_dpsi = d_dpsi(Bcov.theta, psi);
    dBzeta_dpsi  = d_dpsi(Bcov.zeta, psi);
    dBzeta_dth   = d_dtheta(Bcov.zeta, theta);

    Jcon.psi   =  dBzeta_dth ./ Jacobian_abs / mu0;
    Jcon.theta = -dBzeta_dpsi ./ Jacobian_abs / mu0;
    Jcon.zeta  = (dBtheta_dpsi - dBpsi_dth) ./ Jacobian_abs / mu0;

    Jcov.psi   = g11.*Jcon.psi + g12.*Jcon.theta + g13.*Jcon.zeta;
    Jcov.theta = g12.*Jcon.psi + g22.*Jcon.theta + g23.*Jcon.zeta;
    Jcov.zeta  = g13.*Jcon.psi + g23.*Jcon.theta + g33.*Jcon.zeta;

    G = struct();

    G.Jacobian_abs = Jacobian_abs;
    G.Jacobian = Jacobian;

    G.g_cov = struct('psipsi',g11,'psitheta',g12,'psizeta',g13, ...
                     'thetatheta',g22,'thetazeta',g23,'zetazeta',g33);

    G.g_contra = struct( ...
        'psipsi',     gcon11, ...
        'psitheta',   gcon12, ...
        'psizeta',    gcon13, ...
        'thetatheta', gcon22, ...
        'thetazeta',  gcon23, ...
        'zetazeta',   gcon33);

    G.detg = detg;
    
    G.B_contra = Bcon;
    G.B_cov = Bcov;
    G.B2 = B2_direct;
    G.B2_direct = B2_direct;
    G.B2_metric = B2_metric;

    G.J_contra = Jcon;
    G.J_cov = Jcov;
end

function Apsi = d_dpsi(A, psi)
    Apsi = zeros(size(A));

    Apsi(2:end-1,:) = (A(3:end,:) - A(1:end-2,:)) ./ ...
                      (psi(3:end) - psi(1:end-2));

    Apsi(1,:) = (A(2,:) - A(1,:)) ./ (psi(2) - psi(1));
    Apsi(end,:) = (A(end,:) - A(end-1,:)) ./ ...
                  (psi(end) - psi(end-1));
end

function Ath = d_dtheta(A, theta)
    dtheta = theta(2) - theta(1);
    Ath = (circshift(A, -1, 2) - circshift(A, 1, 2)) ./ (2*dtheta);
end