function merluc = CShape_Mercier_Luc(surf, param, bnd)
    R = bnd.R(:);
    Z = bnd.Z(:);
    theta_geo = bnd.theta(:);

    Rt = bnd.Rt(:);
    Zt = bnd.Zt(:);
    Rtt = bnd.Rtt(:);
    Ztt = bnd.Ztt(:);

    dl_dt = hypot(Rt, Zt);

    if any(~isfinite(dl_dt)) || any(dl_dt <= 0)
        error('CShape_Mercier_Luc:DegenerateBoundary', ...
            'The C-shape boundary has a zero or invalid tangent length.');
    end

    cosu_raw =  Rt./dl_dt;
    sinu_raw = -Zt./dl_dt;
    invRc_raw = (Rt.*Ztt - Zt.*Rtt)./dl_dt.^3;

    [Bp, dpsi_dr, Rr, Zr, jac] = CShape_poloidal_field( ...
        theta_geo, R, Rt, Zt, surf.F, surf.q, param);

    % The C-shape parameterization is not required to have the same
    % orientation as the extracted GEQDSK contour.  Choose the Mercier
    % normal sign from its projection on the increasing-r direction.
    normal_projection_raw = sinu_raw.*Rr + cosu_raw.*Zr;
    
    radial_normal_sign = sign(median(normal_projection_raw));
    flux_radial_sign = sign(dpsi_dr);
    
    normal_sign = flux_radial_sign*radial_normal_sign;
    if normal_sign == 0
        error('CShape_Mercier_Luc:UndefinedNormal', ...
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

    % (R,Z) points and coordinate derivatives.
    merluc.R = R;
    merluc.Z = Z;
    merluc.Rt = Rt;
    merluc.Zt = Zt;
    merluc.Rr = Rr;
    merluc.Zr = Zr;

    % Poloidal-coordinate Jacobian and arc-length derivative.
    merluc.jac = jac;
    merluc.dl_dt = dl_dt;

    % Mercier normal and curvature.
    merluc.cosu = cosu;
    merluc.sinu = sinu;
    merluc.cosu_raw = cosu_raw;
    merluc.sinu_raw = sinu_raw;
    merluc.normal_sign = normal_sign;
    merluc.normal_projection_raw = normal_projection_raw;
    merluc.invRc = invRc;
    merluc.invRc_raw = invRc_raw;

    % Magnetic field.
    merluc.dpsi_dr = dpsi_dr;
    merluc.Bp = Bp;
    merluc.Bphi = Bphi;
    merluc.B2 = B2;
    merluc.B = B;
    merluc.H = H;

    % Surface constants and an independent q consistency check.
    merluc.F = surf.F;
    merluc.q = surf.q;
    merluc.q_check = surf.F/(2*pi)*trapz( ...
        theta_geo, dl_dt./(R.^2.*Bp));
end

function [Bp, dpsi_dr, Rr, Zr, jac] = CShape_poloidal_field( ...
        theta, R, Rt, Zt, F, q, param)
    Rr = param.dA_dr ...
        + param.dB_dr*sin(theta) ...
        + param.dC_dr*cos(2*theta);
    Zr = param.dG_dr*cos(theta) ...
        - param.dH_dr*sin(2*theta);

    jac = Rr.*Zt - Rt.*Zr;
    jac_scale = max(1, max(abs(Rr.*Zt)) + max(abs(Rt.*Zr)));
    if any(~isfinite(jac)) || any(abs(jac) <= 1.e-12*jac_scale)
        error('CShape_Mercier_Luc:SingularJacobian', ...
            ['The fitted radial coordinate has a singular Jacobian. ', ...
             'Check the radial coefficient fit and its fitting width.']);
    end

    if ~isfinite(F) || ~isfinite(q) || abs(q) <= eps
        error('CShape_Mercier_Luc:InvalidSurfaceConstant', ...
            'F and q must be finite, and q must be nonzero.');
    end

    q_integral = trapz(theta, abs(jac)./R);
    dpsi_dr_abs = abs(F*q_integral/(2*pi*q));
    
    flux_sign = sign(param.dpsi_dr_fit);
    if flux_sign == 0
        error('CShape_Mercier_Luc:UndefinedFluxDirection', ...
            'Could not determine the sign of dpsi/dr.');
    end
    
    dpsi_dr = flux_sign*dpsi_dr_abs;
    
    Bp = abs(dpsi_dr).*hypot(Rt, Zt)./(abs(jac).*R);
end
