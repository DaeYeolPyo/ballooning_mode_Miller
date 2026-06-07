function bal = eq_ballooning_coefficients(eq, psin0, varargin)
%EQ_BALLOONING_COEFFICIENTS Ballooning-equation coefficients on one flux surface.
%
%   bal = EQ_BALLOONING_COEFFICIENTS(eq, psin0)
%   bal = EQ_BALLOONING_COEFFICIENTS(eq, psin0, 'Name', value, ...)
%
% Requires:
%   - eq_straight_fieldline_theta.m on MATLAB path
%   - eq struct already read from GEQDSK
%
% Main formulas implemented:
%   alpha_t = zeta - q(psin)*(theta - theta0)
%
%   g = (b.grad_N theta) * |grad_N alpha_t|^2 / (B/BN)
%   f = |grad_N alpha_t|^2 / ((b.grad_N theta)*(B/BN)^3)
%   c = 2/(b.grad_N theta)/(B/BN)^4 * d(mu0*p/BN^2)/dpsin ...
%       * [ b x grad_N((2*mu0*p+B^2)/(2*BN^2)) ].grad_N alpha_t
%
% Notes:
%   - psin is used as psi_N.
%   - theta is the straight-field-line theta returned by
%     eq_straight_fieldline_theta with JacMode=4 by default.
%   - zeta is not explicitly needed for axisymmetric coefficients.
%   - For local s-alpha scans, S_hat overrides dq/dpsin and Alpha overrides
%     dp/dpsin through alpha = AlphaFactor*(mu0*dp/dpsin).

    p = inputParser;
    addParameter(p, 'JacMode', 4, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'NTheta', 256, @(x)isnumeric(x) && isscalar(x) && x >= 16);
    addParameter(p, 'Dpsin', 1e-3, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'Theta0', 0.0, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'Alpha0', 0.0, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'S_hat', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'Alpha', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'ShearDefinition', 'flux', @(x)ischar(x) || isstring(x));
    addParameter(p, 'AlphaFactor', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x ~= 0));
    addParameter(p, 'aN', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'BN', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'Mu0', 4*pi*1e-7, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'InterpMethod', 'spline', @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    opt = p.Results;
    opt.InterpMethod = char(lower(string(opt.InterpMethod)));
    if strcmp(opt.InterpMethod, 'pchip')
        opt.InterpMethod = 'spline';   % 2D griddedInterpolant does not support pchip
    end

    if psin0 <= opt.Dpsin || psin0 >= 1-opt.Dpsin
        error('Choose psin0 safely inside the plasma, e.g. Dpsin < psin0 < 1-Dpsin.');
    end

    ntheta = round(opt.NTheta);
    theta = linspace(0, 2*pi, ntheta+1).';
    theta(end) = [];
    dtheta = 2*pi/ntheta;

    psim = psin0 - opt.Dpsin;
    psip = psin0 + opt.Dpsin;

    % ---------------------------------------------------------------
    % 1) Get R(psin,theta), Z(psin,theta) on three nearby surfaces
    % ---------------------------------------------------------------
    cm = surface_on_uniform_theta(eq, psim, theta, opt);
    c0 = surface_on_uniform_theta(eq, psin0, theta, opt);
    cp = surface_on_uniform_theta(eq, psip, theta, opt);

    Rm = cm.R; Zm = cm.Z;
    R  = c0.R; Z  = c0.Z;
    Rp = cp.R; Zp = cp.Z;

    % ---------------------------------------------------------------
    % 2) Flux-coordinate derivatives and gradients in cylindrical basis
    %    basis components are [e_R, e_phi, e_Z]
    % ---------------------------------------------------------------
    R_psin = (Rp - Rm) ./ (psip - psim);
    Z_psin = (Zp - Zm) ./ (psip - psim);

    R_theta = periodic_derivative(R, dtheta);
    Z_theta = periodic_derivative(Z, dtheta);

    D = R_psin .* Z_theta - R_theta .* Z_psin;

    grad_psin  = [ Z_theta ./ D, zeros(ntheta,1), -R_theta ./ D ];
    grad_theta = [-Z_psin ./ D, zeros(ntheta,1),  R_psin ./ D ];
    grad_zeta  = [ zeros(ntheta,1), 1./R, zeros(ntheta,1) ];

    % ---------------------------------------------------------------
    % 3) Profiles: F=R*Bphi, p, q and their psin derivatives
    % ---------------------------------------------------------------
    prof = get_profiles(eq);
    F0  = interp1(prof.psin, prof.F, psin0, 'pchip', 'extrap');
    p0  = interp1(prof.psin, prof.p, psin0, 'pchip', 'extrap');
    q0  = interp1(prof.psin, prof.q, psin0, 'pchip', 'extrap');
    dp0 = interp1(prof.psin, prof.dp_dpsin, psin0, 'pchip', 'extrap');
    dq0 = interp1(prof.psin, prof.dq_dpsin, psin0, 'pchip', 'extrap');
    if ~isempty(opt.S_hat)
        dq0 = shear_to_dq_dpsin(opt.S_hat, q0, psin0, opt.ShearDefinition);
    end

    if isempty(opt.aN)
        if isfield(eq, 'rbbbs')
            aN = 0.5*(max(eq.rbbbs(:)) - min(eq.rbbbs(:)));
        else
            aN = 1.0;
        end
    else
        aN = opt.aN;
    end

    if isempty(opt.BN)
        BN = abs(F0 / eq.rmaxis);
    else
        BN = opt.BN;
    end

    % ---------------------------------------------------------------
    % 4) Magnetic field and B^2 on the central surface
    % ---------------------------------------------------------------
    dpsi_phys = eq.sibry - eq.simag;
    psi_R = dpsi_phys .* grad_psin(:,1);
    psi_Z = dpsi_phys .* grad_psin(:,3);

    BR   = -psi_Z ./ R;
    Bphi =  F0 ./ R;
    BZ   =  psi_R ./ R;

    Bvec = [BR, Bphi, BZ];
    Bmag = sqrt(sum(Bvec.^2, 2));
    bvec = Bvec ./ Bmag;
    BoverBN = Bmag ./ BN;

    alphaFactorDefault = default_alpha_factor(q0, eq, grad_psin, Bmag, BN);
    if isempty(opt.AlphaFactor)
        alphaFactor = alphaFactorDefault;
    else
        alphaFactor = opt.AlphaFactor;
    end
    if ~isempty(opt.Alpha)
        dp0 = (opt.Alpha ./ alphaFactor) ./ opt.Mu0;
    end

    % ---------------------------------------------------------------
    % 5) Clebsch alpha_t gradient
    %    alpha_t = zeta - q(psin)*(theta - theta0)
    % ---------------------------------------------------------------
    theta_shift = theta - opt.Theta0;
    gradN_psin  = aN .* grad_psin;
    gradN_theta = aN .* grad_theta;
    gradN_zeta  = aN .* grad_zeta;
    gradN_alpha = gradN_zeta - q0.*gradN_theta - dq0.*theta_shift.*gradN_psin;
    b_dot_gradN_theta = sum(bvec .* gradN_theta, 2);
    gradN_alpha_sq = sum(gradN_alpha.^2, 2);

    % ---------------------------------------------------------------
    % 6) grad_N H, H=(2*mu0*p+B^2)/(2*BN^2)
    %    Use GEQDSK-grid interpolation for B^2(psin,theta), then finite diff.
    % ---------------------------------------------------------------
    B2m = B2_on_points(eq, Rm, Zm, prof, psim, opt);
    B20 = Bmag.^2;
    B2p = B2_on_points(eq, Rp, Zp, prof, psip, opt);

    B2_psin  = (B2p - B2m) ./ (psip - psim);
    B2_theta = periodic_derivative(B20, dtheta);

    grad_B2 = B2_psin .* grad_psin + B2_theta .* grad_theta;
    grad_H  = (opt.Mu0 * dp0 / BN^2) .* grad_psin + (0.5 / BN^2) .* grad_B2;
    gradN_H = aN .* grad_H;

    dpbar_dpsin = opt.Mu0 * dp0 / BN^2;

    cross_term = sum(cross(bvec, gradN_H, 2) .* gradN_alpha, 2);

    % ---------------------------------------------------------------
    % 7) Ballooning coefficients
    % ---------------------------------------------------------------
    gcoef = b_dot_gradN_theta .* gradN_alpha_sq ./ BoverBN;
    fcoef = gradN_alpha_sq ./ (b_dot_gradN_theta .* BoverBN.^3);
    ccoef = 2 ./ b_dot_gradN_theta ./ BoverBN.^4 .* dpbar_dpsin .* cross_term;

    % ---------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------
    bal = struct();
    bal.psin = psin0;
    bal.theta = theta;
    bal.R = R;
    bal.Z = Z;
    bal.alpha_t_zeta0 = mod(-q0*(theta - opt.Theta0), 2*pi);

    bal.g = gcoef;
    bal.c = ccoef;
    bal.f = fcoef;

    bal.q = q0;
    bal.dq_dpsin = dq0;
    bal.p = p0;
    bal.dp_dpsin = dp0;
    bal.dpbar_dpsin = dpbar_dpsin;
    if isempty(opt.S_hat)
        bal.s_hat = dq0 .* 2 .* psin0 ./ q0;
    else
        bal.s_hat = opt.S_hat;
    end
    bal.alpha = opt.Alpha;
    if isempty(bal.alpha)
        bal.alpha = alphaFactor .* opt.Mu0 .* dp0;
    end
    bal.alphaFactor = alphaFactor;
    bal.alphaFactorDefault = alphaFactorDefault;
    bal.shearDefinition = char(lower(string(opt.ShearDefinition)));
    bal.F = F0;
    bal.BN = BN;
    bal.aN = aN;
    bal.theta0 = opt.Theta0;
    bal.alpha0 = opt.Alpha0;

    bal.BR = BR;
    bal.Bphi = Bphi;
    bal.BZ = BZ;
    bal.bvec = bvec;
    bal.B = Bmag;
    bal.BoverBN = BoverBN;
    bal.b_dot_gradN_theta = b_dot_gradN_theta;
    bal.gradN_alpha_sq = gradN_alpha_sq;
    bal.gradN_alpha = gradN_alpha;
    bal.gradN_psin = gradN_psin;
    bal.gradN_theta = gradN_theta;
    bal.gradN_zeta = gradN_zeta;
    bal.gradN_H = gradN_H;
    bal.cross_term = cross_term;
end

% =====================================================================
function c = surface_on_uniform_theta(eq, psin, theta_grid, opt)
    cont = eq_straight_fieldline_theta(eq, psin, ...
        'JacMode', opt.JacMode, ...
        'NTheta', opt.NTheta, ...
        'InterpMethod', char(opt.InterpMethod));

    th = cont.theta_coord(:);
    R  = cont.R(:);
    Z  = cont.Z(:);

    th = mod(th - th(1), 2*pi);
    [th, ord] = sort(th);
    R = R(ord);
    Z = Z(ord);
    [th, R, Z] = unique_theta(th, R, Z);

    th_ext = [th-2*pi; th; th+2*pi];
    R_ext  = [R; R; R];
    Z_ext  = [Z; Z; Z];

    c.R = interp1(th_ext, R_ext, theta_grid, 'pchip');
    c.Z = interp1(th_ext, Z_ext, theta_grid, 'pchip');
end

% =====================================================================
function dfdtheta = periodic_derivative(f, dtheta)
    f = f(:);
    dfdtheta = (circshift(f,-1) - circshift(f,1)) ./ (2*dtheta);
end

% =====================================================================
function [theta_u, R_u, Z_u] = unique_theta(theta, R, Z)
    theta = theta(:); R = R(:); Z = Z(:);
    valid = isfinite(theta) & isfinite(R) & isfinite(Z);
    theta = theta(valid); R = R(valid); Z = Z(valid);

    if isempty(theta)
        theta_u = theta;
        R_u = R;
        Z_u = Z;
        return;
    end

    tol = 1e-12 * max(1, max(abs(theta)));
    group = cumsum([1; abs(diff(theta)) > tol]);
    n = group(end);
    theta_u = zeros(n, 1);
    R_u = zeros(n, 1);
    Z_u = zeros(n, 1);
    for k = 1:n
        m = group == k;
        theta_u(k) = mean(theta(m));
        R_u(k) = mean(R(m));
        Z_u(k) = mean(Z(m));
    end
end

% =====================================================================
function prof = get_profiles(eq)
    F = get_first_field(eq, {'fpol','ffprime_fpol','F','f','g'});
    p = get_first_field(eq, {'pres','pressure','p'});
    q = get_first_field(eq, {'qpsi','q','qprof'});

    F = F(:); p = p(:); q = q(:);
    n = numel(q);
    psin = linspace(0, 1, n).';

    if numel(F) ~= n
        psinF = linspace(0,1,numel(F)).';
        F = interp1(psinF, F, psin, 'pchip', 'extrap');
    end
    if numel(p) ~= n
        psinp = linspace(0,1,numel(p)).';
        p = interp1(psinp, p, psin, 'pchip', 'extrap');
    end

    prof = struct();
    prof.psin = psin;
    prof.F = F;
    prof.p = p;
    prof.q = q;
    prof.dp_dpsin = gradient(p, psin);
    prof.dq_dpsin = gradient(q, psin);
end

% =====================================================================
function val = get_first_field(s, names)
    for k = 1:numel(names)
        if isfield(s, names{k})
            val = s.(names{k});
            return;
        end
    end
    error('Required profile field not found. Tried: %s', strjoin(names, ', '));
end

% =====================================================================
function dq = shear_to_dq_dpsin(sHat, q0, psin0, definition)
    definition = lower(string(definition));
    switch definition
        case "flux"
            dq = sHat .* q0 ./ max(2 .* psin0, eps);
        case "minor"
            dq = sHat .* q0;
        otherwise
            error('eq_ballooning_coefficients:BadShearDefinition', ...
                'ShearDefinition must be ''flux'' or ''minor''.');
    end
end

% =====================================================================
function alphaFactor = default_alpha_factor(q0, eq, grad_psin, Bmag, BN)
    if isfield(eq, 'rmaxis') && isfinite(eq.rmaxis)
        R0 = eq.rmaxis;
    elseif isfield(eq, 'rcentr') && isfinite(eq.rcentr)
        R0 = eq.rcentr;
    else
        R0 = 1.0;
    end

    gradRef = mean(vecnorm(grad_psin, 2, 2), 'omitnan');
    if ~isfinite(gradRef) || gradRef <= 0
        gradRef = max(vecnorm(grad_psin, 2, 2), [], 'omitnan');
    end

    if ~isfinite(BN) || BN <= 0
        B2ref = mean(Bmag.^2, 'omitnan');
    else
        B2ref = BN.^2;
    end

    alphaFactor = -2 .* q0.^2 .* R0 .* gradRef ./ max(B2ref, eps);
end

% =====================================================================
function B2 = B2_on_points(eq, R, Z, prof, psin, opt)
    [Rgrid, Zgrid, psiRZ] = local_grids_for_B(eq);
    % psiRZ is nR x nZ, so columns are Z and rows are R.
    [psi_Z, psi_R] = gradient(psiRZ, Zgrid, Rgrid);

    psiRF = griddedInterpolant({Rgrid, Zgrid}, psi_R, opt.InterpMethod, 'none');
    psiZF = griddedInterpolant({Rgrid, Zgrid}, psi_Z, opt.InterpMethod, 'none');

    F0 = interp1(prof.psin, prof.F, psin, 'pchip', 'extrap');

    psiR = psiRF(R, Z);
    psiZ = psiZF(R, Z);

    BR = -psiZ ./ R;
    BZ =  psiR ./ R;
    Bphi = F0 ./ R;

    B2 = BR.^2 + BZ.^2 + Bphi.^2;
end

% =====================================================================
function [Rgrid, Zgrid, psiRZ] = local_grids_for_B(eq)
    if isfield(eq, 'rgrid')
        Rgrid = eq.rgrid(:);
    else
        nw = size(eq.psirz, 1);
        if isfield(eq, 'nw'), nw = eq.nw; end
        Rgrid = eq.rleft + (0:nw-1).' * eq.rdim / (nw - 1);
    end

    if isfield(eq, 'zgrid')
        Zgrid = eq.zgrid(:);
    else
        nh = size(eq.psirz, 2);
        if isfield(eq, 'nh'), nh = eq.nh; end
        Zgrid = eq.zmid - eq.zdim/2 + (0:nh-1).' * eq.zdim / (nh - 1);
    end

    psi = eq.psirz;
    if isequal(size(psi), [numel(Rgrid), numel(Zgrid)])
        psiRZ = psi;
    elseif isequal(size(psi), [numel(Zgrid), numel(Rgrid)])
        psiRZ = psi.';
    else
        error('eq.psirz size does not match R/Z grid.');
    end
end
