function [Rt, Zt] = theta_derv(coeffs, ntheta)
    theta = linspace(0, pi*2, ntheta);

    A = coeffs(1, 1);
    B = coeffs(2, 1);
    C = coeffs(3, 1);
    G = coeffs(4, 1);
    H = coeffs(5, 1);

    Rt = A + B*cos(theta) - 2*C*sin(2*theta);
    Zt = -G*sin(theta) - 2*H*cos(2*theta);

    Rt = Rt - A;
end