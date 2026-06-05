function bal = miller_ballooning_coefficients(p, varargin)
%MILLER_BALLOONING_COEFFICIENTS Build Gaur-style ballooning coefficients.
%
%   bal = miller_ballooning_coefficients(p)
%   bal = miller_ballooning_coefficients(p,'Name',value,...)
%
% The input p is the Miller D-shape parameter struct returned by
% DshapeMillerParams. The fields p.s_hat and p.alpha are interpreted as the
% Miller magnetic shear and pressure-gradient parameters. The returned struct
% is compatible with solve_ballooning_eigenvalue.m in ../Gaur.

    ip = inputParser;
    addParameter(ip, 'NTheta', 257, @(x)isnumeric(x)&&isscalar(x)&&x>=32);
    addParameter(ip, 'NGeom', 1201, @(x)isnumeric(x)&&isscalar(x)&&x>=128);
    addParameter(ip, 'Theta0', 0.0, @(x)isnumeric(x)&&isscalar(x));
    addParameter(ip, 'DrFrac', 1e-4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(ip, 'Mu0', 4*pi*1e-7, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(ip, varargin{:});
    opt = ip.Results;

    p = fill_required_defaults(p);
    p.theta0 = opt.Theta0;

    r0 = p.r;
    aN = r0;
    BN = p.B0;
    dr = opt.DrFrac*r0;

    ntheta = round(opt.NTheta);
    theta = linspace(0, 2*pi, ntheta+1).';
    theta(end) = [];

    surf0 = miller_surface_on_pest_theta(p, theta, opt.NGeom);

    [pprime, dpbar_dpsin, geomScales] = miller_pressure_gradient(p, surf0, opt.Mu0);
    FFprime = miller_FFprime_from_shear(p, surf0, pprime, opt.Mu0);

    pm = perturb_miller_surface(p, -dr, FFprime);
    pp = perturb_miller_surface(p,  dr, FFprime);
    surfm = miller_surface_on_pest_theta(pm, theta, opt.NGeom);
    surfp = miller_surface_on_pest_theta(pp, theta, opt.NGeom);

    dpsin = dr/aN;
    R_psin = (surfp.R - surfm.R)./(2*dpsin);
    Z_psin = (surfp.Z - surfm.Z)./(2*dpsin);
    R_theta = periodic_derivative(theta, surf0.R);
    Z_theta = periodic_derivative(theta, surf0.Z);

    D = R_psin.*Z_theta - R_theta.*Z_psin;
    grad_psin  = [ Z_theta./D, zeros(ntheta,1), -R_theta./D ];
    grad_theta = [-Z_psin./D, zeros(ntheta,1),  R_psin./D ];
    grad_zeta  = [zeros(ntheta,1), 1./surf0.R, zeros(ntheta,1)];

    dpsi_dpsin = surf0.dpdr*aN;
    grad_psi = dpsi_dpsin.*grad_psin;
    BR = -grad_psi(:,3)./surf0.R;
    BZ =  grad_psi(:,1)./surf0.R;
    Bphi = surf0.F./surf0.R;
    Bvec = [BR, Bphi, BZ];
    Bmag = sqrt(sum(Bvec.^2, 2));
    bvec = Bvec./Bmag;
    BoverBN = Bmag./BN;

    B2_psin = (surfp.B2 - surfm.B2)./(2*dpsin);
    B20 = Bmag.^2;
    B2_theta = periodic_derivative(theta, B20);
    grad_B2 = B2_psin.*grad_psin + B2_theta.*grad_theta;
    grad_H = dpbar_dpsin.*grad_psin + (0.5/BN^2).*grad_B2;

    gradN_psin = aN.*grad_psin;
    gradN_theta = aN.*grad_theta;
    gradN_zeta = aN.*grad_zeta;
    gradN_H = aN.*grad_H;

    q0 = p.q;
    dq_dpsin = q0*p.s_hat;
    theta_shift = theta - opt.Theta0;
    % Miller varies s and alpha by changing p' and FF' on the local
    % Mercier-Luc equilibrium. In the circular limit this gives the
    % Greene-Chance local shear q*(s*theta - alpha*sin(theta)).
    miller_shear = q0 .* (p.s_hat .* theta_shift - p.alpha .* sin(theta));
    gradN_alpha = gradN_zeta - q0.*gradN_theta - miller_shear.*gradN_psin;
    b_dot_gradN_theta = sum(bvec.*gradN_theta, 2);
    gradN_alpha_sq = sum(gradN_alpha.^2, 2);
    cross_term = sum(cross(bvec, gradN_H, 2).*gradN_alpha, 2);

    gcoef = b_dot_gradN_theta.*gradN_alpha_sq./BoverBN;
    fcoef = gradN_alpha_sq./(b_dot_gradN_theta.*BoverBN.^3);
    ccoef = 2./b_dot_gradN_theta./BoverBN.^4.*dpbar_dpsin.*cross_term;

    bal = struct();
    bal.model = 'miller-local';
    bal.theta = theta;
    bal.R = surf0.R;
    bal.Z = surf0.Z;
    bal.g = gcoef;
    bal.c = ccoef;
    bal.f = fcoef;
    bal.q = q0;
    bal.dq_dpsin = dq_dpsin;
    bal.local_shear_model = 'miller-s-alpha';
    bal.local_shear_periodic_alpha = q0 .* p.alpha;
    bal.s_hat = p.s_hat;
    bal.alpha = p.alpha;
    bal.pprime = pprime;
    bal.FFprime = FFprime;
    bal.dpbar_dpsin = dpbar_dpsin;
    bal.BN = BN;
    bal.aN = aN;
    bal.theta0 = opt.Theta0;
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
    bal.miller = p;
    bal.miller_scales = geomScales;
    bal.surface = surf0;
end

function p = fill_required_defaults(p)
    if ~isfield(p, 'r'), p.r = 1.0; end
    if ~isfield(p, 'B0'), p.B0 = 1.0; end
    if ~isfield(p, 'A'), p.A = 3.0; end
    if ~isfield(p, 'kappa'), p.kappa = 1.0; end
    if ~isfield(p, 'delta'), p.delta = 0.0; end
    if ~isfield(p, 's_kappa'), p.s_kappa = 0.0; end
    if ~isfield(p, 's_delta'), p.s_delta = 0.0; end
    if ~isfield(p, 'dR0_dr'), p.dR0_dr = 0.0; end
    if ~isfield(p, 'q'), p.q = 1.0; end
    if ~isfield(p, 's_hat'), p.s_hat = 0.0; end
    if ~isfield(p, 'alpha'), p.alpha = 0.0; end
    if ~isfield(p, 'F'), p.F = p.A*p.r*p.B0; end
end

function surf = miller_surface_on_pest_theta(p, theta_grid, ngeom)
    [Bp, raw] = Bpol_Dshape(p, ngeom);
    u = raw.u(:);
    R = raw.R(:);
    Z = raw.Z(:);
    Ru = raw.Ru(:);
    Zu = raw.Zu(:);
    jac = raw.jac(:);
    Bp = Bp(:);

    pitch = raw.F.*abs(jac)./(raw.dpdr.*R);
    theta_raw = cumtrapz(u, pitch)./p.q;
    theta_raw = theta_raw - theta_raw(1);
    theta_raw = 2*pi*theta_raw./theta_raw(end);

    [theta_u, keep] = unique_monotone(theta_raw);
    R = R(keep);
    Z = Z(keep);
    Ru = Ru(keep);
    Zu = Zu(keep);
    jac = jac(keep);
    Bp = Bp(keep);

    theta_ext = [theta_u-2*pi; theta_u; theta_u+2*pi];
    R_ext = [R; R; R];
    Z_ext = [Z; Z; Z];
    Bp_ext = [Bp; Bp; Bp];
    Ru_ext = [Ru; Ru; Ru];
    Zu_ext = [Zu; Zu; Zu];
    jac_ext = [jac; jac; jac];

    tq = mod(theta_grid(:), 2*pi);
    surf = struct();
    surf.theta = theta_grid(:);
    surf.R = interp1(theta_ext, R_ext, tq, 'pchip');
    surf.Z = interp1(theta_ext, Z_ext, tq, 'pchip');
    surf.Bp = interp1(theta_ext, Bp_ext, tq, 'pchip');
    surf.Ru = interp1(theta_ext, Ru_ext, tq, 'pchip');
    surf.Zu = interp1(theta_ext, Zu_ext, tq, 'pchip');
    surf.jac = interp1(theta_ext, jac_ext, tq, 'pchip');
    surf.F = raw.F;
    surf.dpdr = raw.dpdr;
    surf.q = p.q;
    surf.Bphi = raw.F./surf.R;
    surf.B2 = surf.Bp.^2 + surf.Bphi.^2;
    surf.raw = raw;
end

function [theta_u, keep] = unique_monotone(theta)
    theta = theta(:);
    keep = [true; diff(theta) > 1e-12];
    theta_u = theta(keep);
    keep_idx = find(keep);
    if numel(theta_u) > 1 && abs(theta_u(end) - theta_u(1) - 2*pi) < 1e-10
        keep_idx(end) = [];
        theta_u(end) = [];
    end
    keep = false(size(theta));
    keep(keep_idx) = true;
end

function [pprime, dpbar_dpsin, scales] = miller_pressure_gradient(p, surf, mu0)
    R0 = p.A*p.r;
    area = polyarea(surf.R, surf.Z);
    centroidR = polygon_centroid_R(surf.R, surf.Z);
    V = 2*pi*area*centroidR;

    raw = surf.raw;
    dVdr = 2*pi*trapz(raw.u(:), raw.R(:).*abs(raw.jac(:)));
    dVdpsi = dVdr/surf.dpdr;

    minor_eff = sqrt(max(V/(2*pi^2*R0), eps));
    pprime = -p.alpha*(2*pi)^2/(2*dVdpsi*minor_eff*mu0);
    dpbar_dpsin = mu0/ p.B0^2 * pprime * surf.dpdr * p.r;

    scales = struct();
    scales.V = V;
    scales.dVdr = dVdr;
    scales.dVdpsi = dVdpsi;
    scales.minor_eff = minor_eff;
end

function cR = polygon_centroid_R(R, Z)
    R = R(:);
    Z = Z(:);
    if R(1) ~= R(end) || Z(1) ~= Z(end)
        R = [R; R(1)];
        Z = [Z; Z(1)];
    end
    crossv = R(1:end-1).*Z(2:end) - R(2:end).*Z(1:end-1);
    A = 0.5*sum(crossv);
    if abs(A) < eps
        cR = mean(R(1:end-1));
    else
        cR = sum((R(1:end-1)+R(2:end)).*crossv)/(6*A);
    end
end

function FFprime = miller_FFprime_from_shear(p, surf, pprime, mu0)
    raw = surf.raw;
    u = raw.u(:);
    R = raw.R(:);
    Ru = raw.Ru(:);
    Zu = raw.Zu(:);
    Bp = Bpol_Dshape(p, numel(u));
    Bp = Bp(:);
    F = raw.F;

    du = u(2) - u(1);
    Ruu = periodic_derivative(u, Ru);
    Zuu = periodic_derivative(u, Zu);
    dl_du = hypot(Ru, Zu);
    invRc = (Ru.*Zuu - Zu.*Ruu)./max(dl_du.^3, eps);
    sinMercier = -Zu./max(dl_du, eps);
    dl = dl_du;

    common = dl./max(R.^3.*Bp.^2, eps);
    I0 = trapz(u, common.*(2*invRc - 2*sinMercier./R));
    Ip = trapz(u, common.*(mu0*pprime.*R./Bp));
    If = trapz(u, common.*(1./(R.*Bp)));

    dq_dpsi = p.q*p.s_hat/(p.r*surf.dpdr);
    denom = p.q/F^2 + F/(2*pi)*If;
    numer = dq_dpsi - F/(2*pi)*(I0 + Ip);
    FFprime = numer/denom;
end

function pp = perturb_miller_surface(p, dr, FFprime)
    pp = p;
    r0 = p.r;
    R0 = p.A*r0;
    x0 = asin(p.delta);
    dpdr0 = [];
    if isfield(p, 'dpdr')
        dpdr0 = p.dpdr;
    end

    pp.r = r0 + dr;
    pp.kappa = p.kappa*(1 + p.s_kappa*dr/r0);
    pp.delta = sin(x0 + p.s_delta*dr/r0);
    pp.A = (R0 + p.dR0_dr*dr)/pp.r;
    pp.q = p.q*(1 + p.s_hat*dr/r0);

    if isempty(dpdr0)
        [~, raw] = Bpol_Dshape(p, 513);
        dpdr0 = raw.dpdr;
    end
    dFdr = FFprime/p.F*dpdr0;
    pp.F = p.F + dFdr*dr;
end

function df = periodic_derivative(x, f)
    x = x(:);
    f = f(:);
    n = numel(f);
    h = mean(diff(x));
    df = zeros(n,1);
    if n < 5
        df = gradient(f, h);
        return;
    end
    
    hasDuplicateEndpoint = abs((x(end) - x(1)) - 2*pi) < 10*h*eps(max(1, abs(x(end)))) ...
        && abs(f(end) - f(1)) < 1e-10 * max(1, max(abs(f)));

    if hasDuplicateEndpoint
        df(2:end-1) = (f(3:end) - f(1:end-2))/(2*h);
        df(1) = (f(2) - f(end-1))/(2*h);
        df(end) = df(1);
    else
        df = (circshift(f, -1) - circshift(f, 1)) ./ (2*h);
    end
end
