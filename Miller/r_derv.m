function [Rr, Zr] = r_derv(coeffs, ntheta)
    theta = linspace(0, pi*2, ntheta);

    Ar = coeffs(1, 2);
    Br = coeffs(2, 2);
    Cr = coeffs(3, 2);
    Gr = coeffs(4, 2);
    Hr = coeffs(5, 2);

    Rr = Ar + Br*sin(theta) + Cr*cos(2*theta);
    Zr = Gr*cos(theta) - Hr*sin(2*theta);
end