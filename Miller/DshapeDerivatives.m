function [Ru, Zu, Rr, Zr, jac, u, R, Z, p] = DshapeDerivatives(p, ntheta)
    [R, Z, u, p] = DshapeParam(p, ntheta);

    r = p.r;
    kappa = p.kappa;
    x = asin(p.delta);
    eta = u + x*sin(u);

    dx_dr = p.s_delta/r;
    dkappa_dr = p.s_kappa*kappa/r;

    Ru = -r*sin(eta).*(1 + x*cos(u));
    Zu = kappa*r*cos(u);

    Rr = p.dR0_dr + cos(eta) - r*sin(eta).*dx_dr.*sin(u);
    Zr = (kappa + r*dkappa_dr).*sin(u);

    jac = Rr.*Zu - Ru.*Zr;
end
