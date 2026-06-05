function [Bp, out] = Bpol_Dshape(p, ntheta)
    [Ru, Zu, Rr, Zr, jac, u, R, Z, p] = DshapeDerivatives(p, ntheta);

    R0 = p.A*p.r;
    if isfield(p, 'F')
        F = p.F;
    else
        F = R0*p.B0;
    end

    if isfield(p, 'dpdr')
        dpdr = p.dpdr;
    else
        q_integral = trapz(u, abs(jac)./R);
        dpdr = F*q_integral/(2*pi*p.q);
    end

    Bp = abs(dpdr).*hypot(Ru, Zu)./(abs(jac).*R);

    out = struct();
    out.u = u;
    out.R = R;
    out.Z = Z;
    out.Ru = Ru;
    out.Zu = Zu;
    out.Rr = Rr;
    out.Zr = Zr;
    out.jac = jac;
    out.F = F;
    out.dpdr = dpdr;
    out.q_check = F*trapz(u, abs(jac)./R)/(2*pi*dpdr);
    out.params = p;
end
