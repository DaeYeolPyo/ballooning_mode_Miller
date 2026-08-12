function bnd = Miller_boundary(R0, r, delta, kappa, opts)
    arguments
        R0 (1,1) double
        r (1,1) double
        delta (1,1) double
        kappa (1, 1) double
        opts.NTheta (1,1) int32 = 256
    end

    x = asin(delta);
    theta = linspace(0, 2*pi, opts.NTheta);

    t1 = theta + x*sin(theta);
    t2 = 1 + x*cos(theta);

    R = R0 + r*cos(t1);
    Rt = -r*t2.*sin(t1);
    Rtt = r*x*sin(theta).*sin(t1) - r*(t2.^2).*cos(t1);

    Z = kappa*r*sin(theta);
    Zt = kappa*r*cos(theta);
    Ztt = -kappa*r*sin(theta);

    bnd = struct();
    bnd.R = R;
    bnd.Rt = Rt;
    bnd.Rtt = Rtt;
    bnd.Z = Z;
    bnd.Zt = Zt;
    bnd.Ztt = Ztt;
    bnd.theta = theta;
end