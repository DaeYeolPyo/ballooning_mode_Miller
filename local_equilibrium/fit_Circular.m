function [param, bnd] = fit_Circular(eq, eqfunc, psiN, opts)
    arguments
        eq (1,1) struct
        eqfunc (1,1) struct
        psiN (1,1) double
        opts.NTheta (1,1) int32 = 256
        opts.dpsi (1,1) double = 3.e-2
        opts.NRadialFit (1,1) int32 = 9
    end

    surf = extract_flux_surface(eq, eqfunc, psiN, Npoints=opts.NTheta);
    [R0, r, fit_rms, fit_max] = fit_one_circle(surf);

    nradial = double(opts.NRadialFit);
    if nradial < 3 || mod(nradial, 2) == 0
        error('fit_Circular:BadNRadialFit', ...
            'NRadialFit must be an odd integer greater than or equal to 3.');
    end

    half_width = min(opts.dpsi, 0.95*min(psiN, 1 - psiN));
    if half_width <= 0
        error('fit_Circular:BadPsiN', ...
            'psiN must lie strictly between 0 and 1.');
    end

    psiN_fit = linspace(psiN - half_width, psiN + half_width, nradial).';
    R0_fit = zeros(nradial, 1);
    r_fit = zeros(nradial, 1);
    residual_rms_fit = zeros(nradial, 1);
    residual_max_fit = zeros(nradial, 1);

    center_index = (nradial + 1)/2;
    for k = 1:nradial
        if k == center_index
            fit_surf = surf;
        else
            fit_surf = extract_flux_surface( ...
                eq, eqfunc, psiN_fit(k), Npoints=opts.NTheta);
        end

        [R0_fit(k), r_fit(k), residual_rms_fit(k), ...
            residual_max_fit(k)] = fit_one_circle(fit_surf);
    end

    % The only radial shaping quantity required by the shifted-circle
    % coordinate is dR0/dr.  Fitting R0 directly against the fitted minor
    % radius avoids differentiating the normalized flux label.
    radial_offset = r_fit - r;
    design = [ones(nradial, 1), radial_offset];
    coeff_R0 = design\R0_fit;
    dR0_dr = coeff_R0(2);

    param = struct();
    param.R0 = R0;
    param.r = r;
    param.A = R0/r;
    param.dR0_dr = dR0_dr;
    param.F = interp1(linspace(0, 1, eq.nw), eq.fpol(:), psiN);
    param.q = interp1(linspace(0, 1, eq.nw), eq.qpsi(:), psiN);
    param.B0 = param.F/param.R0;
    param.fit_residual_rms = fit_rms;
    param.fit_residual_max = fit_max;
    param.radial_fit = struct( ...
        'psiN', psiN_fit, ...
        'R0', R0_fit, ...
        'r', r_fit, ...
        'residual_rms', residual_rms_fit, ...
        'residual_max', residual_max_fit, ...
        'half_width', half_width, ...
        'npoints', nradial);

    bnd = Circular_boundary(R0, r, NTheta=opts.NTheta);
    bnd.Req = surf.R;
    bnd.Zeq = surf.Z;
end

function [R0, r, residual_rms, residual_max] = fit_one_circle(surf)
    R = surf.R(:);
    Z = surf.Z(:);

    R_Z0 = find_R_at_Z0_from_surface(surf);
    R0_midplane = mean(R_Z0);
    r_midplane = abs(diff(R_Z0))/2;

    % Minimize the spread of the distance from a center constrained to Z=0.
    % This uses the whole contour and is therefore insensitive to the noisy
    % top-point triangularity that motivated the circular model.
    objective = @(center) mean(( ...
        hypot(R - center, Z) - mean(hypot(R - center, Z))).^2);
    lower = R0_midplane - 0.5*r_midplane;
    upper = R0_midplane + 0.5*r_midplane;
    options = optimset('Display', 'off', 'TolX', 1.e-12);
    R0 = fminbnd(objective, lower, upper, options);

    radial_distance = hypot(R - R0, Z);
    r = mean(radial_distance);
    residual = radial_distance - r;
    residual_rms = sqrt(mean(residual.^2));
    residual_max = max(abs(residual));
end
