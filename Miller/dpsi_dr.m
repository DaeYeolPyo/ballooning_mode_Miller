function dpdr = dpsi_dr(eq, psiN, ntheta)
    psiN_interp = linspace(0, 1, eq.nw);
    fpol = eq.fpol;
    qpsi = abs(eq.qpsi);

    F = spline(psiN_interp, fpol, psiN);
    q = spline(psiN_interp, qpsi, psiN);

    coeffs = Cshape_metrics(eq, psiN);
    [R, ~] = CshapeParam(coeffs(:, 1), ntheta);
    [Rr, Zr] = r_derv(coeffs, ntheta);
    [Rt, Zt] = theta_derv(coeffs, ntheta);

    theta = linspace(0, 2*pi, ntheta);
    integrand = abs(Rr.*Zt - Rt.*Zr)./R;

    int = trapz(theta, integrand);
    dpdr = (F/q)*int/(2*pi);
end
