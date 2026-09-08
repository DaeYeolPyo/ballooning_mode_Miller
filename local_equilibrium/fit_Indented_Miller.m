function [param, bnd] = fit_Indented_Miller(eq, eqfunc, psiN, opts)
    arguments
        eq (1,1) struct
        eqfunc (1,1) struct
        psiN (1,1) double
        opts.NTheta (1,1) int32 = 256
        opts.dpsi (1,1) double = 3.e-2
        opts.NRadialFit (1,1) int32 = 9
        opts.MaxIterations (1,1) int32 = 200
        opts.ThetaTolerance (1,1) double = 1.e-9
        opts.CoefficientTolerance (1,1) double = 1.e-10
        opts.OptimizerTolerance (1,1) double = 1.e-10
        opts.NPhaseStarts (1,1) int32 = 32
        opts.NXStarts (1,1) int32 = 33
        opts.MaxAbsDelta (1,1) double = 0.98
        opts.MaxKappa (1,1) double = 5.0
        opts.MaxI0 (1,1) double = 3.0
        opts.Verbose (1,1) logical = false
    end

    validate_options(psiN, opts);
    nradial = double(opts.NRadialFit);
    half_width = min(opts.dpsi, 0.95*min(psiN, 1 - psiN));
    if ~isfinite(half_width) || half_width <= 0
        error('fit_Indented_Miller:BadRadialFitWidth', ...
            'dpsi must define a finite positive radial fitting interval.');
    end

    surf = extract_flux_surface( ...
        eq, eqfunc, psiN, Npoints=opts.NTheta);
    [coeff, center_fit] = fit_single_surface(surf, opts);

    % The stable fitting coefficients are
    %   [R0, a, x, b, d] = [R0, a, asin(delta), a*kappa, a*I0].
    % Fit their radial variation against a itself.  This avoids the
    % singular-looking ratios that arise if kappa and I0 are differentiated
    % before the geometric length scales b and d.
    psiN_fit = linspace( ...
        psiN - half_width, psiN + half_width, nradial).';
    coeff_fit = zeros(nradial, 5);
    fit_diagnostics = cell(nradial, 1);
    center_index = (nradial + 1)/2;

    for k = 1:nradial
        if k == center_index
            coeff_fit(k, :) = coeff.';
            fit_diagnostics{k} = center_fit;
        else
            fit_surf = extract_flux_surface( ...
                eq, eqfunc, psiN_fit(k), Npoints=opts.NTheta);
            [coeff_k, fit_diagnostics{k}] = ...
                fit_single_surface(fit_surf, opts);
            coeff_fit(k, :) = coeff_k.';
        end
    end

    a = coeff(2);
    radial_offset = coeff_fit(:, 2) - a;
    if max(radial_offset) - min(radial_offset) <= ...
            1.e-10*max(1, abs(a))
        error('fit_Indented_Miller:DegenerateRadialCoordinate', ...
            'The fitted minor radius does not vary across radial surfaces.');
    end

    design = [ones(nradial, 1), radial_offset];
    radial_regression = design\coeff_fit;
    dcoeff_da = radial_regression(2, :).';

    psi_fit = eq.simag + psiN_fit*(eq.sibry - eq.simag);
    radial_flux_regression = design\psi_fit;
    dpsi_da_fit = radial_flux_regression(2);

    R0 = coeff(1);
    x = coeff(3);
    b = coeff(4);
    d = coeff(5);
    delta = sin(x);
    kappa = b/a;
    I0 = d/a;

    param = struct();
    param.R0 = R0;
    param.a = a;
    param.r = a;
    param.A = R0/a;
    param.x = x;
    param.delta = delta;
    param.b = b;
    param.kappa = kappa;
    param.d = d;
    param.I0 = I0;
    param.indentation = I0;
    param.coeff = coeff;

    param.dR0_da = dcoeff_da(1);
    param.dR0_dr = param.dR0_da;
    param.da_da = dcoeff_da(2);
    param.dx_da = dcoeff_da(3);
    param.db_da = dcoeff_da(4);
    param.dd_da = dcoeff_da(5);
    param.dcoeff_da = dcoeff_da;
    param.coeff_table = [coeff, dcoeff_da];

    % Miller-style dimensionless radial shaping derivatives.
    param.s_delta = a*param.dx_da;
    param.s_kappa = a*param.db_da/b - 1;
    param.dI0_da = (param.dd_da - I0)/a;
    param.dI0_dr = param.dI0_da;
    param.s_I0 = a*param.dI0_da;
    param.s_I = param.s_I0;
    param.dpsi_da_fit = dpsi_da_fit;
    param.dpsi_dr_fit = dpsi_da_fit;

    param.F = interp1(linspace(0, 1, eq.nw), eq.fpol(:), psiN);
    param.q = interp1(linspace(0, 1, eq.nw), eq.qpsi(:), psiN);
    param.B0 = param.F/R0;
    param.fit = center_fit;
    param.radial_fit = struct( ...
        'psiN', psiN_fit, ...
        'a', coeff_fit(:, 2), ...
        'R0', coeff_fit(:, 1), ...
        'x', coeff_fit(:, 3), ...
        'delta', sin(coeff_fit(:, 3)), ...
        'b', coeff_fit(:, 4), ...
        'kappa', coeff_fit(:, 4)./coeff_fit(:, 2), ...
        'd', coeff_fit(:, 5), ...
        'I0', coeff_fit(:, 5)./coeff_fit(:, 2), ...
        'coeff', coeff_fit, ...
        'dcoeff_da', dcoeff_da, ...
        'half_width', half_width, ...
        'npoints', nradial, ...
        'diagnostics', {fit_diagnostics});

    bnd = Indented_Miller_boundary( ...
        R0, a, delta, kappa, I0, NTheta=opts.NTheta);
    bnd.Req = surf.R;
    bnd.Zeq = surf.Z;
end

function validate_options(psiN, opts)
    if psiN <= 0 || psiN >= 1
        error('fit_Indented_Miller:BadPsiN', ...
            'psiN must lie strictly between 0 and 1.');
    end
    if opts.NTheta < 16 || opts.MaxIterations < 1 ...
            || opts.NPhaseStarts < 8 || opts.NXStarts < 9
        error('fit_Indented_Miller:BadResolution', ...
            ['NTheta must be at least 16, MaxIterations must be positive, ', ...
             'NPhaseStarts at least 8, and NXStarts at least 9.']);
    end
    if opts.NRadialFit < 3 || mod(double(opts.NRadialFit), 2) == 0
        error('fit_Indented_Miller:BadNRadialFit', ...
            'NRadialFit must be an odd integer at least 3.');
    end
    if opts.MaxAbsDelta <= 0 || opts.MaxAbsDelta >= 1 ...
            || opts.MaxKappa <= 0 || opts.MaxI0 <= 0
        error('fit_Indented_Miller:BadBounds', ...
            ['Require 0 < MaxAbsDelta < 1, MaxKappa > 0, ', ...
             'and MaxI0 > 0.']);
    end
end

function [coeff, diag] = fit_single_surface(surf, opts)
    [R, Z] = prepare_contour(surf.R, surf.Z);
    npoint = numel(R);
    if npoint < 8
        error('fit_Indented_Miller:TooFewSurfacePoints', ...
            'At least eight distinct contour points are required.');
    end

    optimizer_options = optimset( ...
        'Display', 'off', 'TolX', opts.OptimizerTolerance);
    orientations = [1, -1];
    best_objective = inf;
    best = struct();

    % The indentation term fixes theta=0 at the outboard indentation, so a
    % robust phase search is necessary.  Arc length supplies only a smooth
    % ordering, not the analytic poloidal angle itself.
    for orientation = orientations
        if orientation > 0
            R_try = R;
            Z_try = Z;
        else
            R_try = flipud(R);
            Z_try = flipud(Z);
        end

        theta_arc = arc_length_angle(R_try, Z_try);
        phase_step = 2*pi/double(opts.NPhaseStarts);
        phase_grid = (0:double(opts.NPhaseStarts)-1)*phase_step;
        phase_objective = zeros(size(phase_grid));

        for j = 1:numel(phase_grid)
            phase_objective(j) = fixed_phase_objective( ...
                phase_grid(j), theta_arc, R_try, Z_try, opts);
        end

        [~, jbest] = min(phase_objective);
        phase_seed = phase_grid(jbest);
        phase = fminbnd(@(value) fixed_phase_objective( ...
            value, theta_arc, R_try, Z_try, opts), ...
            phase_seed - phase_step, phase_seed + phase_step, ...
            optimizer_options);

        theta_try = theta_arc + phase;
        [coeff_try, objective_try] = fit_coefficients( ...
            theta_try, R_try, Z_try, opts);

        if objective_try < best_objective && is_valid_coeff(coeff_try, opts)
            best_objective = objective_try;
            best.R = R_try;
            best.Z = Z_try;
            best.theta = theta_try;
            best.coeff = coeff_try;
            best.orientation = orientation;
            best.phase = phase;
        end
    end

    if isempty(fieldnames(best))
        error('fit_Indented_Miller:InitializationFailed', ...
            ['No physically admissible phase/orientation initialization ', ...
             'was found. Increase MaxI0 or MaxKappa if appropriate.']);
    end

    R = best.R;
    Z = best.Z;
    theta = best.theta;
    coeff = best.coeff;
    initial_coeff = coeff;
    initial_objective = best_objective;
    objective_history = nan(double(opts.MaxIterations) + 1, 1);
    objective_history(1) = initial_objective;
    converged = false;
    theta_change = inf;
    coefficient_change = inf;

    for iteration = 1:double(opts.MaxIterations)
        theta_new = closest_point_angles(theta, R, Z, coeff);
        [coeff_new, objective_new] = fit_coefficients( ...
            theta_new, R, Z, opts, coeff(3));

        if ~is_valid_coeff(coeff_new, opts)
            break;
        end

        objective_history(iteration + 1) = objective_new;
        theta_change = max(abs(theta_new - theta));
        coefficient_change = norm(coeff_new - coeff) ...
            /max(norm(coeff), eps);
        theta = theta_new;
        coeff = coeff_new;

        if opts.Verbose
            fprintf(['  Indented-Miller iteration %3d: RMS %.6e, ', ...
                'dtheta %.3e, dcoeff %.3e\n'], ...
                iteration, sqrt(objective_new), ...
                theta_change, coefficient_change);
        end

        if theta_change <= opts.ThetaTolerance ...
                && coefficient_change <= opts.CoefficientTolerance
            converged = true;
            break;
        end
    end

    finite_history = isfinite(objective_history);
    objective_history = objective_history(finite_history);
    [Rfit, Zfit] = evaluate_indented_geometry(coeff, theta);
    distance = hypot(R - Rfit, Z - Zfit);
    characteristic_radius = 0.25*((max(R) - min(R)) ...
        + (max(Z) - min(Z)));

    diag = struct();
    diag.converged = converged;
    diag.iterations = iteration;
    diag.rms_distance = sqrt(mean(distance.^2));
    diag.max_distance = max(distance);
    diag.relative_rms_distance = diag.rms_distance ...
        /max(characteristic_radius, eps);
    diag.initial_rms_distance = sqrt(initial_objective);
    diag.objective_history = objective_history;
    diag.final_theta_change = theta_change;
    diag.final_coefficient_change = coefficient_change;
    if numel(objective_history) >= 2
        diag.final_relative_objective_change = ...
            abs(objective_history(end) - objective_history(end - 1)) ...
            /max(objective_history(end - 1), eps);
    else
        diag.final_relative_objective_change = NaN;
    end
    diag.theta = mod(theta, 2*pi);
    diag.orientation = best.orientation;
    diag.initial_phase = mod(best.phase, 2*pi);
    diag.initial_coeff = initial_coeff;
    diag.vertical_center_offset = 0.5*(max(Z) + min(Z));

    if ~converged && opts.Verbose
        warning('fit_Indented_Miller:NoConvergence', ...
            ['Closest-point refinement reached iteration %d. ', ...
             'Final relative RMS error is %.6e.'], ...
            iteration, diag.relative_rms_distance);
    end
end

function theta_new = closest_point_angles(theta_old, R, Z, coeff)
    theta_old = theta_old(:);
    R = R(:);
    Z = Z(:);
    npoint = numel(theta_old);
    period = 2*pi;
    nsearch = max(256, 4*ceil(sqrt(npoint)));
    theta_search = (0:nsearch - 1)*(period/nsearch);
    search_step = period/nsearch;

    [Rsearch, Zsearch] = evaluate_indented_geometry( ...
        coeff, theta_search);
    distance_search = (R - Rsearch).^2 + (Z - Zsearch).^2;
    [~, nearest_index] = min(distance_search, [], 2);
    theta_seed = theta_search(nearest_index).';
    theta = theta_seed + period*round((theta_old - theta_seed)/period);

    max_newton_iterations = 15;
    newton_tolerance = 1.e-11;
    max_newton_step = 2*search_step;

    for iteration = 1:max_newton_iterations
        [Rfit, Zfit, Rt, Zt, Rtt, Ztt] = ...
            evaluate_indented_geometry(coeff, theta);
        dR = Rfit - R;
        dZ = Zfit - Z;
        numerator = dR.*Rt + dZ.*Zt;
        denominator = Rt.^2 + Zt.^2 + dR.*Rtt + dZ.*Ztt;
        denominator_scale = Rt.^2 + Zt.^2 ...
            + abs(dR.*Rtt) + abs(dZ.*Ztt);

        valid = abs(denominator) > ...
            1.e-12*max(denominator_scale, eps);
        step = zeros(npoint, 1);
        step(valid) = numerator(valid)./denominator(valid);
        step = max(-max_newton_step, min(max_newton_step, step));

        old_distance = dR.^2 + dZ.^2;
        step_scale = ones(npoint, 1);
        theta_candidate = theta;
        best_distance = old_distance;

        for backtrack = 1:8
            theta_trial = theta - step_scale.*step;
            [Rtrial, Ztrial] = evaluate_indented_geometry( ...
                coeff, theta_trial);
            trial_distance = (Rtrial - R).^2 + (Ztrial - Z).^2;
            improved = trial_distance < best_distance;
            theta_candidate(improved) = theta_trial(improved);
            best_distance(improved) = trial_distance(improved);
            step_scale(~improved) = 0.5*step_scale(~improved);
        end

        theta_change = max(abs(theta_candidate - theta));
        theta = theta_candidate;
        if theta_change <= newton_tolerance
            break;
        end
    end

    minimum_separation = 64*eps(max(1, max(abs(theta))));
    theta_ordered = project_to_ordered_angles(theta, minimum_separation);
    old_objective = geometric_objective(theta_old, R, Z, coeff);
    theta_new = theta_old;

    for backtrack = 0:10
        fraction = 2^(-backtrack);
        candidate = theta_old + fraction*(theta_ordered - theta_old);
        order_is_valid = all(diff(candidate) > 0) ...
            && candidate(end) < candidate(1) + period;
        if order_is_valid
            candidate_objective = geometric_objective( ...
                candidate, R, Z, coeff);
            if candidate_objective <= old_objective*(1 + 1.e-12)
                theta_new = candidate;
                return;
            end
        end
    end
end

function theta_ordered = project_to_ordered_angles(theta, separation)
    n = numel(theta);
    ramp = (0:n-1).'*separation;
    y = theta(:) - ramp;
    level = zeros(n, 1);
    weight = zeros(n, 1);
    first = zeros(n, 1);
    last = zeros(n, 1);
    nblock = 0;

    for i = 1:n
        nblock = nblock + 1;
        level(nblock) = y(i);
        weight(nblock) = 1;
        first(nblock) = i;
        last(nblock) = i;
        while nblock > 1 && level(nblock - 1) > level(nblock)
            combined_weight = weight(nblock - 1) + weight(nblock);
            level(nblock - 1) = (weight(nblock - 1)*level(nblock - 1) ...
                + weight(nblock)*level(nblock))/combined_weight;
            weight(nblock - 1) = combined_weight;
            last(nblock - 1) = last(nblock);
            nblock = nblock - 1;
        end
    end

    theta_ordered = zeros(n, 1);
    for block = 1:nblock
        theta_ordered(first(block):last(block)) = level(block);
    end
    theta_ordered = theta_ordered + ramp;
end

function value = fixed_phase_objective(phase, theta_arc, R, Z, opts)
    theta = theta_arc + phase;
    [~, value] = fit_coefficients(theta, R, Z, opts);
end

function [coeff, value] = fit_coefficients(theta, R, Z, opts, x_seed)
    xmin = -asin(opts.MaxAbsDelta);
    xmax = asin(opts.MaxAbsDelta);
    nx = double(opts.NXStarts);
    xgrid = linspace(xmin, xmax, nx);

    if nargin < 5
        values = zeros(size(xgrid));
        for j = 1:nx
            [~, values(j)] = coefficients_at_x( ...
                xgrid(j), theta, R, Z, opts);
        end
        [~, jbest] = min(values);
        x_grid_best = xgrid(jbest);
        lower_index = max(1, jbest - 1);
        upper_index = min(nx, jbest + 1);
        lower = xgrid(lower_index);
        upper = xgrid(upper_index);
    else
        x_grid_best = min(xmax, max(xmin, x_seed));
        local_half_width = 2*(xmax - xmin)/(nx - 1);
        lower = max(xmin, x_grid_best - local_half_width);
        upper = min(xmax, x_grid_best + local_half_width);
    end

    optimizer_options = optimset( ...
        'Display', 'off', 'TolX', opts.OptimizerTolerance);
    if upper - lower > opts.OptimizerTolerance
        x_opt = fminbnd(@(x) scalar_objective( ...
            x, theta, R, Z, opts), lower, upper, optimizer_options);
        [coeff_opt, value_opt] = coefficients_at_x( ...
            x_opt, theta, R, Z, opts);
    else
        coeff_opt = nan(5, 1);
        value_opt = inf;
    end

    [coeff_grid, value_grid] = coefficients_at_x( ...
        x_grid_best, theta, R, Z, opts);
    if value_opt < value_grid
        coeff = coeff_opt;
        value = value_opt;
    else
        coeff = coeff_grid;
        value = value_grid;
    end
end

function value = scalar_objective(x, theta, R, Z, opts)
    [~, value] = coefficients_at_x(x, theta, R, Z, opts);
end

function [coeff, value] = coefficients_at_x(x, theta, R, Z, opts)
    theta = theta(:);
    R = R(:);
    Z = Z(:);
    eta = theta + x*sin(theta);
    h = ((1 + cos(theta))/2).^3;

    radial_coeff = [ones(size(theta)), cos(eta), -h]\R;
    b = sin(theta)\Z;
    coeff = [radial_coeff(1); radial_coeff(2); x; b; radial_coeff(3)];

    [Rfit, Zfit] = evaluate_indented_geometry(coeff, theta);
    base_value = mean((R - Rfit).^2 + (Z - Zfit).^2);
    if is_valid_coeff(coeff, opts)
        value = base_value;
    else
        scale = max(1, 0.25*((max(R) - min(R)) ...
            + (max(Z) - min(Z))));
        value = base_value + scale^2*(100 + constraint_penalty(coeff, opts));
    end
end

function penalty = constraint_penalty(coeff, opts)
    R0 = coeff(1);
    a = coeff(2);
    b = coeff(4);
    d = coeff(5);
    safe_a = max(abs(a), eps);
    kappa = b/safe_a;
    I0 = d/safe_a;
    penalty = max(0, -R0/safe_a)^2 ...
        + max(0, -a/safe_a)^2 ...
        + max(0, -b/safe_a)^2 ...
        + max(0, -d/safe_a)^2 ...
        + max(0, kappa - opts.MaxKappa)^2 ...
        + max(0, I0 - opts.MaxI0)^2;
end

function valid = is_valid_coeff(coeff, opts)
    valid = numel(coeff) == 5 && all(isfinite(coeff)) ...
        && coeff(1) > 0 && coeff(2) > 0 && coeff(4) > 0 ...
        && coeff(5) >= 0 ...
        && abs(sin(coeff(3))) <= opts.MaxAbsDelta*(1 + 1.e-12) ...
        && coeff(4)/coeff(2) <= opts.MaxKappa ...
        && coeff(5)/coeff(2) <= opts.MaxI0;
end

function value = geometric_objective(theta, R, Z, coeff)
    [Rfit, Zfit] = evaluate_indented_geometry(coeff, theta);
    value = mean((R(:) - Rfit(:)).^2 + (Z(:) - Zfit(:)).^2);
end

function [R, Z, Rt, Zt, Rtt, Ztt] = ...
        evaluate_indented_geometry(coeff, theta)
    R0 = coeff(1);
    a = coeff(2);
    x = coeff(3);
    b = coeff(4);
    d = coeff(5);

    cosine = cos(theta);
    sine = sin(theta);
    eta = theta + x*sine;
    eta_t = 1 + x*cosine;
    eta_tt = -x*sine;
    h = ((1 + cosine)/2).^3;

    R = R0 + a*cos(eta) - d*h;
    Z = b*sine;

    if nargout >= 3
        h_t = -(3/8)*(1 + cosine).^2.*sine;
        Rt = -a*sin(eta).*eta_t - d*h_t;
        Zt = b*cosine;
    end

    if nargout >= 5
        h_tt = (3/4)*(1 + cosine).*sine.^2 ...
            - (3/8)*(1 + cosine).^2.*cosine;
        Rtt = -a*cos(eta).*eta_t.^2 ...
            - a*sin(eta).*eta_tt - d*h_tt;
        Ztt = -b*sine;
    end
end

function theta = arc_length_angle(R, Z)
    Rclosed = [R(:); R(1)];
    Zclosed = [Z(:); Z(1)];
    ds = hypot(diff(Rclosed), diff(Zclosed));
    perimeter = sum(ds);
    if ~isfinite(perimeter) || perimeter <= 0
        error('fit_Indented_Miller:DegenerateContour', ...
            'The input contour has zero or non-finite perimeter.');
    end
    s = [0; cumsum(ds(1:end-1))];
    theta = 2*pi*s/perimeter;
end

function [R, Z] = prepare_contour(R, Z)
    R = R(:);
    Z = Z(:);
    if numel(R) ~= numel(Z) || any(~isfinite([R; Z]))
        error('fit_Indented_Miller:BadContour', ...
            'Surface R and Z must be finite vectors of equal length.');
    end

    scale = max(1, max(hypot(R, Z)));
    if numel(R) > 1 && hypot(R(end) - R(1), Z(end) - Z(1)) ...
            <= 1.e-12*scale
        R(end) = [];
        Z(end) = [];
    end

    step = hypot(circshift(R, -1) - R, circshift(Z, -1) - Z);
    if any(step <= 16*eps(scale))
        error('fit_Indented_Miller:DuplicateContourPoint', ...
            'The contour contains coincident neighboring points.');
    end
end
