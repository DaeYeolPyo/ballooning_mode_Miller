function bnd = Indented_Miller_boundary(R0, a, delta, kappa, I0, opts)
    arguments
        R0 (1,1) double
        a (1,1) double
        delta (1,1) double
        kappa (1,1) double
        I0 (1,1) double
        opts.NTheta (1,1) int32 = 256
    end

    if ~isfinite(R0) || ~isfinite(a) || ~isfinite(delta) ...
            || ~isfinite(kappa) || ~isfinite(I0)
        error('Indented_Miller_boundary:NonFiniteInput', ...
            'All shape parameters must be finite.');
    end
    if a <= 0 || kappa <= 0 || abs(delta) >= 1 || I0 < 0
        error('Indented_Miller_boundary:BadShape', ...
            ['Require a > 0, kappa > 0, |delta| < 1, ', ...
             'and I0 >= 0.']);
    end
    if opts.NTheta < 16
        error('Indented_Miller_boundary:BadResolution', ...
            'NTheta must be at least 16.');
    end

    theta = linspace(0, 2*pi, double(opts.NTheta));
    x = asin(delta);
    b = kappa*a;
    d = I0*a;

    cosine = cos(theta);
    sine = sin(theta);
    eta = theta + x*sine;
    eta_t = 1 + x*cosine;
    eta_tt = -x*sine;

    h = ((1 + cosine)/2).^3;
    h_t = -(3/8)*(1 + cosine).^2.*sine;
    h_tt = (3/4)*(1 + cosine).*sine.^2 ...
        - (3/8)*(1 + cosine).^2.*cosine;

    R = R0 + a*cos(eta) - d*h;
    Rt = -a*sin(eta).*eta_t - d*h_t;
    Rtt = -a*cos(eta).*eta_t.^2 ...
        - a*sin(eta).*eta_tt - d*h_tt;

    Z = b*sine;
    Zt = b*cosine;
    Ztt = -b*sine;

    bnd = struct();
    bnd.R = R;
    bnd.Rt = Rt;
    bnd.Rtt = Rtt;
    bnd.Z = Z;
    bnd.Zt = Zt;
    bnd.Ztt = Ztt;
    bnd.theta = theta;
    bnd.eta = eta;
    bnd.h = h;
    bnd.x = x;
    bnd.b = b;
    bnd.d = d;
end
