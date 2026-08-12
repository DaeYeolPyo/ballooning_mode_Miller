function merluc = Miller_Mercier_Luc(surf, param, bnd)
    R = bnd.R(:);
    Z = bnd.Z(:);
    theta_geo = bnd.theta(:);

    Rt = bnd.Rt(:);
    Zt = bnd.Zt(:);
    Rtt = bnd.Rtt(:);
    Ztt = bnd.Ztt(:);

    dl_dt = hypot(Rt, Zt);

    cosu_raw =  Rt./dl_dt;
    sinu_raw = -Zt./dl_dt;
    invRc_raw = (Rt.*Ztt - Zt.*Rtt)./dl_dt.^3;

    [Bp, dpsi_dr, Rr, Zr, jac] = Miller_poloidal_field( ...
        theta_geo, R, Rt, Zt, surf.F, surf.q, param);

    % Align the Mercier normal with increasing Miller r (and increasing psi
    % for the present equilibrium).  The CCW boundary parameterization makes
    % (sinu_raw, cosu_raw) point inward, so its sign must not be used blindly
    % in the radial equilibrium identities.
    normal_projection_raw = sinu_raw.*Rr + cosu_raw.*Zr;
    normal_sign = sign(median(normal_projection_raw));
    if normal_sign == 0
        error('Miller_Mercier_Luc:UndefinedNormal', ...
            'Could not determine the Mercier normal orientation.');
    end

    cosu = normal_sign*cosu_raw;
    sinu = normal_sign*sinu_raw;
    invRc = normal_sign*invRc_raw;

    Bphi = surf.F./R;
    B2 = Bp.^2 + Bphi.^2;
    B = sqrt(B2);

    H = surf.F./(surf.q*(R.^2).*Bp);

    theta_PEST = cumtrapz(theta_geo, H.*dl_dt);

    merluc = struct();
    % poloidal angles
    merluc.theta_geo = theta_geo;
    merluc.theta_PEST = theta_PEST;

    % (R, Z) Cartesian points/derivatives
    merluc.R = R;
    merluc.Z = Z;
    merluc.Rt = Rt;
    merluc.Zt = Zt;
    merluc.Rr = Rr;
    merluc.Zr = Zr;

    % Jacobian
    merluc.jac = jac;

    % Poloidal arc
    merluc.dl_dt = dl_dt;

    % u angle
    merluc.cosu = cosu;
    merluc.sinu = sinu;
    merluc.cosu_raw = cosu_raw;
    merluc.sinu_raw = sinu_raw;
    merluc.normal_sign = normal_sign;
    merluc.normal_projection_raw = normal_projection_raw;

    % Curvature radius
    merluc.invRc = invRc;
    merluc.invRc_raw = invRc_raw;

    % dpsi/dr & poloidal field
    merluc.dpsi_dr = dpsi_dr;
    merluc.Bp = Bp;

    % toroidal field & total magnetic field strength
    merluc.Bphi = Bphi;
    merluc.B2 = B2;
    merluc.B = B;

    %
    merluc.H = H;

    % surface constants
    merluc.F = surf.F;
    merluc.q = surf.q;
    merluc.q_check = surf.F/(2*pi)*trapz( ...
        theta_geo, dl_dt./(R.^2.*Bp));
end

function [Bp, dpsi_dr, Rr, Zr, jac] = Miller_poloidal_field( ...
        theta, R, Rt, Zt, F, q, param)
    delta = param.delta;
    kappa = param.kappa;
    r = param.r;

    x = asin(delta);
    eta = theta + x*sin(theta);

    dx_dr = param.s_delta/r;
    dkappa_dr = param.s_kappa*kappa/r;

    Rr = param.dR0_dr + cos(eta) ...
        - r*sin(eta).*dx_dr.*sin(theta);
    Zr = (kappa + r*dkappa_dr).*sin(theta);

    jac = Rr.*Zt - Rt.*Zr;

    q_integral = trapz(theta, abs(jac)./R);
    dpsi_dr = F*q_integral/(2*pi*q);

    Bp = abs(dpsi_dr).*hypot(Rt, Zt)./(abs(jac).*R);
end
