function circluc = Circular_Mercier_Luc(surf, param, bnd)
    R = bnd.R(:);
    Z = bnd.Z(:);
    theta_geo = bnd.theta(:);

    Rt = bnd.Rt(:);
    Zt = bnd.Zt(:);
    Rtt = bnd.Rtt(:);
    Ztt = bnd.Ztt(:);
    dl_dt = hypot(Rt, Zt);

    cosu_raw = Rt./dl_dt;
    sinu_raw = -Zt./dl_dt;
    invRc_raw = (Rt.*Ztt - Zt.*Rtt)./dl_dt.^3;

    Rr = param.dR0_dr + cos(theta_geo);
    Zr = sin(theta_geo);
    radial_metric_factor = 1 + param.dR0_dr*cos(theta_geo);
    if any(radial_metric_factor <= 0)
        error('Circular_Mercier_Luc:SingularRadialCoordinate', ...
            ['The shifted-circle coordinate requires ', ...
             '1 + dR0/dr*cos(theta) > 0 everywhere.']);
    end

    jac = Rr.*Zt - Rt.*Zr;
    jac_closed_form = param.r.*radial_metric_factor;
    jac_relative_error = max(abs(jac - jac_closed_form)) ...
        /max(max(abs(jac_closed_form)), eps);

    % Exact shifted-circle q integral:
    % integral |J|/R dtheta = 2*pi*[r/D + R0'*(1 - R0/D)],
    % where D=sqrt(R0^2-r^2).  The coordinate regularity check above makes
    % the absolute value unnecessary.
    if param.R0 <= param.r
        error('Circular_Mercier_Luc:InvalidAspectRatio', ...
            'The circular toroidal surface requires R0 > r.');
    end
    D = sqrt(param.R0^2 - param.r^2);
    q_integral_analytic = 2*pi*( ...
        param.r/D + param.dR0_dr*(1 - param.R0/D));
    q_integral_numeric = trapz(theta_geo, jac_closed_form./R);
    q_integral_relative_error = abs( ...
        q_integral_numeric - q_integral_analytic) ...
        /max(abs(q_integral_analytic), eps);

    dpsi_dr = surf.F*q_integral_analytic/(2*pi*surf.q);
    Bp = abs(dpsi_dr)./(R.*radial_metric_factor);

    normal_projection_raw = sinu_raw.*Rr + cosu_raw.*Zr;
    normal_sign = sign(median(normal_projection_raw));
    if normal_sign == 0
        error('Circular_Mercier_Luc:UndefinedNormal', ...
            'Could not determine the circular-surface normal orientation.');
    end

    cosu = normal_sign*cosu_raw;
    sinu = normal_sign*sinu_raw;
    invRc = normal_sign*invRc_raw;

    Bphi = surf.F./R;
    B2 = Bp.^2 + Bphi.^2;
    B = sqrt(B2);
    H = surf.F./(surf.q*(R.^2).*Bp);
    theta_PEST = cumtrapz(theta_geo, H.*dl_dt);

    circluc = struct();
    circluc.theta_geo = theta_geo;
    circluc.theta_PEST = theta_PEST;
    circluc.R = R;
    circluc.Z = Z;
    circluc.Rt = Rt;
    circluc.Zt = Zt;
    circluc.Rr = Rr;
    circluc.Zr = Zr;
    circluc.jac = jac_closed_form;
    circluc.dl_dt = dl_dt;
    circluc.cosu = cosu;
    circluc.sinu = sinu;
    circluc.cosu_raw = cosu_raw;
    circluc.sinu_raw = sinu_raw;
    circluc.normal_sign = normal_sign;
    circluc.normal_projection_raw = normal_projection_raw;
    circluc.invRc = invRc;
    circluc.invRc_raw = invRc_raw;
    circluc.dpsi_dr = dpsi_dr;
    circluc.Bp = Bp;
    circluc.Bphi = Bphi;
    circluc.B2 = B2;
    circluc.B = B;
    circluc.H = H;
    circluc.F = surf.F;
    circluc.q = surf.q;
    circluc.q_check = surf.F/(2*pi)*trapz( ...
        theta_geo, dl_dt./(R.^2.*Bp));
    circluc.radial_metric_factor = radial_metric_factor;
    circluc.q_integral_analytic = q_integral_analytic;
    circluc.q_integral_numeric = q_integral_numeric;
    circluc.q_integral_relative_error = q_integral_relative_error;
    circluc.jac_relative_error = jac_relative_error;
end
