function [R, Z] = CshapeParam(coeffs, ntheta)
    theta = linspace(0, pi*2, ntheta);

    A = coeffs(1);
    B = coeffs(2);
    C = coeffs(3);
    G = coeffs(4);
    H = coeffs(5);

    R = A + B*sin(theta) + C*cos(2*theta);
    Z = G*cos(theta) - H*sin(2*theta);
end