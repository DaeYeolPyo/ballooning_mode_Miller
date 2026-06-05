function param = fit_Miller(eq, psiN, varargin)
    p = inputParser;
    addParameter(p, 'NTheta', 256, @(x)isnumeric(x) && isscalar(x) && x >= 8);
    addParameter(p, 'dpsi', 1.e-6, @(x)isfloat && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    dp = opt.dpsi;
    surf = extract_flux_surface(eq, psiN, 'NTheta', opt.NTheta);
    surfp = extract_flux_surface(eq, psiN+dp, 'NTheta', opt.NTheta);
    surfm = extract_flux_surface(eq, psiN-dp, 'NTheta', opt.NTheta);

    [R0, r, delta, kappa] = find_shaping(surf);
    [R01, r1, delta1, kappa1] = find_shaping(surfp);
    [R02, r2, delta2, kappa2] = find_shaping(surfm);

    drR0 = (R02 - R01)/(r2 - r1);
    drk = (kappa2 - kappa1)/(r2 - r1);
    drd = (delta2 - delta1)/(r2 - r1);

    s_kappa = r*drk/kappa;
    s_delta = r*drd/sqrt(1 - delta^2);

    param = struct();
    param.R0 = R0;
    param.r = r;
    param.A = R0/r;
    param.delta = delta;
    param.kappa = kappa;
    param.drR0 = drR0;
    param.s_delta = s_delta;
    param.s_kappa = s_kappa;
end

function [R0, r, delta, kappa] = find_shaping(s)
    R_Z0 = find_R_at_Z0_from_surface(s);

    R0 = (R_Z0(2) + R_Z0(1))*0.5;
    r = abs(R_Z0(2) - R_Z0(1))*0.5;
    
    [Zmax, thZInd] = max(s.Z);
    Rtop = s.R(thZInd);

    delta = -(Rtop - R0)/r;
    kappa = Zmax/r;
end