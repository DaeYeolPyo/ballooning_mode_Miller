function [param, bnd] = fit_Miller(eq, eqfunc, psiN, opts)
    arguments
        eq (1,1) struct
        eqfunc (1,1) struct
        psiN (1,1) double
        opts.NTheta (1,1) int32 = 256
        opts.dpsi (1,1) double = 3.e-2
        opts.NRadialFit (1,1) int32 = 9
    end

    surf = extract_flux_surface(eq, eqfunc, psiN, Npoints=opts.NTheta);
    [R0, r, delta, kappa] = find_shaping(surf, opts.NTheta);

    nradial = double(opts.NRadialFit);
    if nradial < 3 || mod(nradial, 2) == 0
        error('fit_Miller:BadNRadialFit', ...
            'NRadialFit must be an odd integer greater than or equal to 3.');
    end

    half_width = min(opts.dpsi, 0.95*min(psiN, 1 - psiN));
    if half_width <= 0
        error('fit_Miller:BadPsiN', ...
            'psiN must lie strictly between 0 and 1.');
    end

    psiN_fit = linspace(psiN - half_width, psiN + half_width, nradial).';
    R0_fit = zeros(nradial, 1);
    r_fit = zeros(nradial, 1);
    delta_fit = zeros(nradial, 1);
    kappa_fit = zeros(nradial, 1);

    center_index = (nradial + 1)/2;
    for k = 1:nradial
        if k == center_index
            fit_surf = surf;
        else
            fit_surf = extract_flux_surface( ...
                eq, eqfunc, psiN_fit(k), Npoints=opts.NTheta);
        end

        [R0_fit(k), r_fit(k), delta_fit(k), kappa_fit(k)] = ...
            find_shaping(fit_surf, opts.NTheta);
    end

    % A local linear regression is deliberately used here.  Differencing two
    % nearly coincident contours amplifies the GEQDSK grid-scale jitter in the
    % top-point triangularity.  Fit x = asin(delta), since s_delta = r*dx/dr.
    radial_offset = r_fit - r;
    design = [ones(nradial, 1), radial_offset];
    coeff_R0 = design\R0_fit;
    coeff_kappa = design\kappa_fit;
    coeff_x = design\asin(max(-1, min(1, delta_fit)));

    drR0 = coeff_R0(2);
    drk = coeff_kappa(2);
    dx_dr = coeff_x(2);

    s_kappa = r*drk/kappa;
    s_delta = r*dx_dr;

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
    param.radial_fit = struct( ...
        'psiN', psiN_fit, ...
        'R0', R0_fit, ...
        'r', r_fit, ...
        'delta', delta_fit, ...
        'kappa', kappa_fit, ...
        'half_width', half_width, ...
        'npoints', nradial);

    bnd = Miller_boundary(R0, r, delta, kappa, Ntheta=opts.NTheta);
    bnd.Req = surf.R;
    bnd.Zeq = surf.Z;
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