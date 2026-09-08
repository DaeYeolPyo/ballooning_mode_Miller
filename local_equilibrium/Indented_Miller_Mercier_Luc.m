function merluc = Indented_Miller_Mercier_Luc(surf, param, bnd)
    R = bnd.R(:);
    Z = bnd.Z(:);
    theta_geo = bnd.theta(:);

    Rt = bnd.Rt(:);
    Zt = bnd.Zt(:);
    Rtt = bnd.Rtt(:);
    Ztt = bnd.Ztt(:);
    dl_dt = hypot(Rt, Zt);

    if any(~isfinite(dl_dt)) || any(dl_dt <= 0)
        error('Indented_Miller_Mercier_Luc:DegenerateBoundary', ...
            ['The indented-Miller boundary has a zero or invalid ', ...
             'tangent length.']);
    end

    cosu_raw = Rt./dl_dt;
    sinu_raw = -Zt./dl_dt;
    invRc_raw = (Rt.*Ztt - Zt.*Rtt)./dl_dt.^3;

    [Bp, dpsi_da, Ra, Za, jac] = indented_poloidal_field( ...
        theta_geo, R, Rt, Zt, surf.F, surf.q, param);

    % The raw Mercier normal is fixed by the direction of increasing theta.
    % Align the final normal with grad(psi), allowing either sign of dpsi/da.
    normal_projection_raw = sinu_raw.*Ra + cosu_raw.*Za;
    radial_normal_sign = sign(median(normal_projection_raw));
    flux_radial_sign = sign(dpsi_da);
    normal_sign = flux_radial_sign*radial_normal_sign;
    if normal_sign == 0
        error('Indented_Miller_Mercier_Luc:UndefinedNormal', ...
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

    % Poloidal angles.
    merluc.theta_geo = theta_geo;
    merluc.theta_PEST = theta_PEST;

    % Analytic surface and local radial-coordinate derivatives.
    merluc.R = R;
    merluc.Z = Z;
    merluc.Rt = Rt;
    merluc.Zt = Zt;
    merluc.Ra = Ra;
    merluc.Za = Za;
    merluc.Rr = Ra;
    merluc.Zr = Za;

    % Poloidal-coordinate Jacobian and arc length.
    merluc.jac = jac;
    merluc.dl_dt = dl_dt;

    % Mercier normal and curvature, aligned with grad(psi).
    merluc.cosu = cosu;
    merluc.sinu = sinu;
    merluc.cosu_raw = cosu_raw;
    merluc.sinu_raw = sinu_raw;
    merluc.normal_sign = normal_sign;
    merluc.normal_projection_raw = normal_projection_raw;
    merluc.invRc = invRc;
    merluc.invRc_raw = invRc_raw;

    % Magnetic field.
    merluc.dpsi_da = dpsi_da;
    merluc.dpsi_dr = dpsi_da;
    merluc.Bp = Bp;
    merluc.Bphi = Bphi;
    merluc.B2 = B2;
    merluc.B = B;
    merluc.H = H;

    % Surface constants and independent safety-factor reconstruction.
    merluc.F = surf.F;
    merluc.q = surf.q;
    merluc.q_check = surf.F/(2*pi)*trapz( ...
        theta_geo, dl_dt./(R.^2.*Bp));
end

function [Bp, dpsi_da, Ra, Za, jac] = indented_poloidal_field( ...
        theta, R, Rt, Zt, F, q, param)
    eta = theta + param.x*sin(theta);
    h = ((1 + cos(theta))/2).^3;

    % Partial derivatives at fixed analytic poloidal angle theta.  The
    % fitted coefficient vector is [R0, a, x, b, d], with a as the radial
    % coordinate and d=a*I0, b=a*kappa.
    Ra = param.dR0_da ...
        + param.da_da*cos(eta) ...
        - param.a*sin(eta).*param.dx_da.*sin(theta) ...
        - param.dd_da*h;
    Za = param.db_da*sin(theta);

    jac = Ra.*Zt - Rt.*Za;
    jac_scale = max(1, max(abs(Ra.*Zt)) + max(abs(Rt.*Za)));
    if any(~isfinite(jac)) || any(abs(jac) <= 1.e-12*jac_scale)
        error('Indented_Miller_Mercier_Luc:SingularJacobian', ...
            ['The fitted (a,theta) coordinate has a singular Jacobian. ', ...
             'Check the radial fit interval and neighboring contours.']);
    end

    if ~isfinite(F) || ~isfinite(q) || abs(q) <= eps
        error('Indented_Miller_Mercier_Luc:InvalidSurfaceConstant', ...
            'F and q must be finite, and q must be nonzero.');
    end

    q_integral = trapz(theta, abs(jac)./R);
    dpsi_da_abs = abs(F*q_integral/(2*pi*q));
    flux_sign = sign(param.dpsi_da_fit);
    if flux_sign == 0
        error('Indented_Miller_Mercier_Luc:UndefinedFluxDirection', ...
            'Could not determine the sign of dpsi/da.');
    end
    dpsi_da = flux_sign*dpsi_da_abs;

    Bp = abs(dpsi_da).*hypot(Rt, Zt)./(abs(jac).*R);
end
