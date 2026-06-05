function Bp = Bpol_from_coeffs(coeffs, dpdr, ntheta)
    [R, ~] = CshapeParam(coeffs(:, 1), ntheta);
    [Rr, Zr] = r_derv(coeffs, ntheta);
    [Rt, Zt] = theta_derv(coeffs, ntheta);

    jac = Rr.*Zt - Rt.*Zr;
    Bp = abs(dpdr).*hypot(Rt, Zt)./(abs(jac).*R);
end
