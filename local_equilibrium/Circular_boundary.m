function bnd = Circular_boundary(R0, r, opts)
    arguments
        R0 (1,1) double
        r (1,1) double {mustBePositive}
        opts.NTheta (1,1) int32 = 256
    end

    theta = linspace(0, 2*pi, opts.NTheta);

    R = R0 + r*cos(theta);
    Rt = -r*sin(theta);
    Rtt = -r*cos(theta);

    Z = r*sin(theta);
    Zt = r*cos(theta);
    Ztt = -r*sin(theta);

    bnd = struct();
    bnd.theta = theta;
    bnd.R = R;
    bnd.Z = Z;
    bnd.Rt = Rt;
    bnd.Zt = Zt;
    bnd.Rtt = Rtt;
    bnd.Ztt = Ztt;
end
