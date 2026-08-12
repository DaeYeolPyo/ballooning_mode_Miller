function bnd = CShape_boundary(A, B, C, G, H, opts)
    arguments
        A (1,1) double
        B (1,1) double
        C (1,1) double
        G (1,1) double
        H (1,1) double
        opts.NTheta (1,1) int32 = 256
    end

    theta = linspace(0, 2*pi, opts.NTheta);

    R = A + B*sin(theta) + C*cos(2*theta);
    Rt = B*cos(theta) - 2*C*sin(2*theta);
    Rtt = -B*sin(theta) - 4*C*cos(2*theta);

    Z = G*cos(theta) - H*sin(2*theta);
    Zt = -G*sin(theta) - 2*H*cos(2*theta);
    Ztt = -G*cos(theta) + 4*H*sin(2*theta);

    bnd = struct();
    bnd.R = R;
    bnd.Rt = Rt;
    bnd.Rtt = Rtt;
    bnd.Z = Z;
    bnd.Zt = Zt;
    bnd.Ztt = Ztt;
    bnd.theta = theta;
end
