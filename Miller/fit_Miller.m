function param = fit_Miller(eq, psiN, varargin)
    p = inputParser;
    addParameter(p, 'NTheta', 256, @(x)isnumeric(x) && isscalar(x) && x >= 8);
    addParameter(p, 'dpsi', 1.e-3, @(x)isfloat(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    dp = opt.dpsi;
    surf = extract_flux_surface(eq, psiN, 'NTheta', opt.NTheta);
    surfp = extract_flux_surface(eq, psiN+dp, 'NTheta', opt.NTheta);
    surfm = extract_flux_surface(eq, psiN-dp, 'NTheta', opt.NTheta);

    [R0, r, delta, kappa] = find_shaping(surf, opt.NTheta);
    [R01, r1, delta1, kappa1] = find_shaping(surfp, opt.NTheta);
    [R02, r2, delta2, kappa2] = find_shaping(surfm, opt.NTheta);

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
    param.dR0_dr = drR0;
    param.s_delta = s_delta;
    param.s_kappa = s_kappa;
    param.F = interp1(linspace(0, 1, eq.nw), eq.fpol(:), psiN);
    param.q = interp1(linspace(0, 1, eq.nw), eq.qpsi(:), psiN);
    param.B0 = param.F/param.R0;
end

function [R0, r, delta, kappa] = find_shaping(s, ntheta)
    R_Z0 = find_R_at_Z0_from_surface(s);

    R0 = (R_Z0(2) + R_Z0(1))*0.5;
    r = abs(R_Z0(2) - R_Z0(1))*0.5;

    [Rtop, Zmax] = find_top_point(s, ntheta);

    delta = -(Rtop - R0)/r;
    kappa = Zmax/r;
end

function [Rtop, Ztop] = find_top_point(s, ntheta)
    R = s.R(:);
    Z = s.Z(:);

    if hypot(R(end) - R(1), Z(end) - Z(1)) > 1.e-12
        R = [R; R(1)];
        Z = [Z; Z(1)];
    end

    arc = cumsum([0; hypot(diff(R), diff(Z))]);
    [arc, idx] = unique(arc, 'stable');
    R = R(idx);
    Z = Z(idx);

    if numel(arc) < 4 || arc(end) == arc(1)
        [Ztop, ind] = max(Z);
        Rtop = R(ind);
        return;
    end

    nsample = max(round(4*ntheta), 2048);
    aq = linspace(arc(1), arc(end), nsample).';
    Zq = interp1(arc, Z, aq, 'spline');
    [~, ind] = max(Zq);

    if ind == 1 || ind == nsample
        atop = aq(ind);
    else
        zfun = @(a) -interp1(arc, Z, a, 'spline');
        atop = fminbnd(zfun, aq(ind-1), aq(ind+1));
    end

    Rtop = interp1(arc, R, atop, 'spline');
    Ztop = interp1(arc, Z, atop, 'spline');
end
