function [param, bnd] = fit_CShape(eq, eqfunc, psiN, opts)
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
        opts.Verbose (1,1) logical = false
    end

    if psiN <= 0 || psiN >= 1
        error('fit_CShape:BadPsiN', ...
            'psiN must lie strictly between 0 and 1.');
    end
    if opts.NTheta < 16 || opts.MaxIterations < 1 ...
            || opts.NPhaseStarts < 8
        error('fit_CShape:BadResolution', ...
            ['NTheta must be at least 16, MaxIterations must be positive, ', ...
             'and NPhaseStarts must be at least 8.']);
    end

    nradial = double(opts.NRadialFit);
    if nradial < 3 || mod(nradial, 2) == 0
        error('fit_CShape:BadNRadialFit', ...
            'NRadialFit must be an odd integer greater than or equal to 3.');
    end

    half_width = min(opts.dpsi, 0.95*min(psiN, 1 - psiN));
    if ~isfinite(half_width) || half_width <= 0
        error('fit_CShape:BadRadialFitWidth', ...
            'dpsi must define a finite positive radial fitting interval.');
    end

    surf = extract_flux_surface( ...
        eq, eqfunc, psiN, Npoints=opts.NTheta);
    [coeff, center_fit] = fit_single_surface(surf, opts);

    % Fit the radial variation of the five C-shape coefficients.  B is the
    % natural midplane half-width of this parameterization and is used as
    % the local minor-radius coordinate r.
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

    r = coeff(2);
    radial_offset = coeff_fit(:, 2) - r;
    if max(radial_offset) - min(radial_offset) <= ...
            1.e-10*max(1, abs(r))
        error('fit_CShape:DegenerateRadialCoordinate', ...
            'The fitted B coefficient does not vary across radial surfaces.');
    end

    design = [ones(nradial, 1), radial_offset];
    radial_regression = design\coeff_fit;
    dcoeff_dr = radial_regression(2, :).';

    psi_fit = eq.simag ...
    + psiN_fit*(eq.sibry - eq.simag);

    radial_flux_regression = design\psi_fit;
    dpsi_dr_fit = radial_flux_regression(2);

    A = coeff(1);
    B = coeff(2);
    C = coeff(3);
    G = coeff(4);
    H = coeff(5);
    R0 = A - C;

    param = struct();
    param.A = A;
    param.B = B;
    param.C = C;
    param.G = G;
    param.H = H;
    param.coeff = coeff;
    param.r = r;
    param.R0 = R0;
    param.aspect_ratio = R0/r;

    param.dA_dr = dcoeff_dr(1);
    param.dB_dr = dcoeff_dr(2);
    param.dC_dr = dcoeff_dr(3);
    param.dG_dr = dcoeff_dr(4);
    param.dH_dr = dcoeff_dr(5);
    param.dcoeff_dr = dcoeff_dr;
    param.coeff_table = [coeff, dcoeff_dr];

    param.F = interp1(linspace(0, 1, eq.nw), eq.fpol(:), psiN);
    param.q = interp1(linspace(0, 1, eq.nw), eq.qpsi(:), psiN);
    param.B0 = param.F/R0;
    param.fit = center_fit;
    param.dpsi_dr_fit = dpsi_dr_fit;
    param.radial_fit = struct( ...
        'psiN', psiN_fit, ...
        'r', coeff_fit(:, 2), ...
        'coeff', coeff_fit, ...
        'dcoeff_dr', dcoeff_dr, ...
        'half_width', half_width, ...
        'npoints', nradial, ...
        'diagnostics', {fit_diagnostics});

    bnd = CShape_boundary(A, B, C, G, H, NTheta=opts.NTheta);
    bnd.Req = surf.R;
    bnd.Zeq = surf.Z;
end

function [coeff, diag] = fit_single_surface(surf, opts)
    [R, Z] = prepare_contour(surf.R, surf.Z);
    npoint = numel(R);

    if npoint < 8
        error('fit_CShape:TooFewSurfacePoints', ...
            'At least eight distinct contour points are required.');
    end

    optimizer_options = optimset( ...
        'Display', 'off', 'TolX', opts.OptimizerTolerance);

    % The restricted Fourier basis fixes both phase and orientation.  Find a
    % robust arc-length initialization by testing the two contour directions
    % and optimizing the common phase shift for each one.
    orientations = [1, -1];
    best_objective = inf;
    best = struct();

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
                phase_grid(j), theta_arc, R_try, Z_try);
        end

        [~, jbest] = min(phase_objective);
        phase_seed = phase_grid(jbest);
        phase = fminbnd(@(value) fixed_phase_objective( ...
            value, theta_arc, R_try, Z_try), ...
            phase_seed - phase_step, phase_seed + phase_step, ...
            optimizer_options);

        theta_try = theta_arc + phase;
        coeff_try = linear_coefficients(theta_try, R_try, Z_try);
        objective_try = geometric_objective( ...
            theta_try, R_try, Z_try, coeff_try);

        if objective_try < best_objective
            best_objective = objective_try;
            best.R = R_try;
            best.Z = Z_try;
            best.theta = theta_try;
            best.coeff = coeff_try;
            best.orientation = orientation;
            best.phase = phase;
        end
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

    % Alternating minimization:
    %   1. coefficients fixed -> closest-point theta_i (argmin),
    %   2. theta_i fixed      -> linear least-squares coefficients.
    for iteration = 1:double(opts.MaxIterations)
        theta_new = closest_point_angles(theta, R, Z, coeff);
        coeff_new = linear_coefficients(theta_new, R, Z);
        objective_new = geometric_objective( ...
            theta_new, R, Z, coeff_new);
        objective_history(iteration + 1) = objective_new;

        theta_change = max(abs(theta_new - theta));
        coefficient_change = norm(coeff_new - coeff) ...
            /max(norm(coeff), eps);

        theta = theta_new;
        coeff = coeff_new;

        if opts.Verbose
            fprintf(['  C-shape iteration %2d: RMS %.6e, ', ...
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

    objective_history = objective_history(1:iteration + 1);
    [coeff, canonical_transform] = canonicalize_coefficients(coeff);

    % Evaluate the final geometric error with the uncanonicalized theta;
    % canonicalization only re-labels the same closed curve.
    coeff_for_correspondence = linear_coefficients(theta, R, Z);
    [Rfit, Zfit] = evaluate_cshape(coeff_for_correspondence, theta);
    distance = hypot(R - Rfit, Z - Zfit);
    characteristic_radius = 0.25*(range(R) + range(Z));

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
    diag.canonical_transform = canonical_transform;
    diag.vertical_center_offset = 0.5*(max(Z) + min(Z));

    if ~converged && opts.Verbose
        warning('fit_CShape:NoConvergence', ...
            ['Closest-point refinement reached MaxIterations=%d. ', ...
             'Final relative RMS error is %.6e.'], ...
            opts.MaxIterations, diag.relative_rms_distance);
    end
end

function theta_new = closest_point_angles(theta_old, R, Z, coeff)
    % Find all point-to-curve projections simultaneously.  A coarse global
    % search selects the correct basin on a possibly concave C-shape, then a
    % damped vectorized Newton iteration solves
    %   (x(theta)-x_i) dot x_theta(theta) = 0.
    theta_old = theta_old(:);
    R = R(:);
    Z = Z(:);
    npoint = numel(theta_old);
    period = 2*pi;
    nsearch = max(128, 2*ceil(sqrt(npoint)));
    theta_search = (0:nsearch - 1)*(period/nsearch);
    search_step = period/nsearch;

    [Rsearch, Zsearch] = evaluate_cshape_geometry( ...
        coeff, theta_search);
    distance_search = (R - Rsearch).^2 + (Z - Zsearch).^2;
    [~, nearest_index] = min(distance_search, [], 2);
    theta_seed = theta_search(nearest_index).';

    % Put every periodic seed on the same unwrapped branch as its previous
    % correspondence.  This makes the ordering projection below meaningful.
    theta = theta_seed + period*round( ...
        (theta_old - theta_seed)/period);

    max_newton_iterations = 12;
    newton_tolerance = 1.e-11;
    max_newton_step = 2*search_step;

    for iteration = 1:max_newton_iterations
        [Rfit, Zfit, Rt, Zt, Rtt, Ztt] = ...
            evaluate_cshape_geometry(coeff, theta);

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

        % Per-point backtracking remains vectorized.  A point accepts a
        % Newton step only when its geometric distance decreases.
        for backtrack = 1:8
            theta_trial = theta - step_scale.*step;
            [Rtrial, Ztrial] = evaluate_cshape_geometry( ...
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

    % Reconcile independently projected points with the monotone phase
    % constraint.  If that projection raises the total objective, damp the
    % update toward the previously valid correspondence.
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
    % Pool-adjacent-violators projection onto a strictly increasing sequence.
    % Subtracting a tiny ramp turns strict ordering into ordinary isotonic
    % ordering; the ramp is restored after the projection.
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

function value = fixed_phase_objective(phase, theta_arc, R, Z)
    theta = theta_arc + phase;
    coeff = linear_coefficients(theta, R, Z);
    value = geometric_objective(theta, R, Z, coeff);
end

function value = geometric_objective(theta, R, Z, coeff)
    [Rfit, Zfit] = evaluate_cshape(coeff, theta);
    value = mean((R - Rfit).^2 + (Z - Zfit).^2);
end

function coeff = linear_coefficients(theta, R, Z)
    theta = theta(:);
    XR = [ones(size(theta)), sin(theta), cos(2*theta)];
    XZ = [cos(theta), -sin(2*theta)];

    coeff_R = XR\R(:);
    coeff_Z = XZ\Z(:);
    coeff = [coeff_R; coeff_Z];
end

function [R, Z] = evaluate_cshape(coeff, theta)
    [R, Z] = evaluate_cshape_geometry(coeff, theta);
end

function [R, Z, Rt, Zt, Rtt, Ztt] = ...
        evaluate_cshape_geometry(coeff, theta)
    A = coeff(1);
    B = coeff(2);
    C = coeff(3);
    G = coeff(4);
    H = coeff(5);

    R = A + B*sin(theta) + C*cos(2*theta);
    Z = G*cos(theta) - H*sin(2*theta);

    if nargout >= 3
        Rt = B*cos(theta) - 2*C*sin(2*theta);
        Zt = -G*sin(theta) - 2*H*cos(2*theta);
    end

    if nargout >= 5
        Rtt = -B*sin(theta) - 4*C*cos(2*theta);
        Ztt = -G*cos(theta) + 4*H*sin(2*theta);
    end
end

function [coeff, transform] = canonicalize_coefficients(coeff)
    % Four phase/orientation symmetries describe the same geometric curve.
    % Choose the conventional representation B >= 0 and G >= 0 so radial
    % fits do not acquire arbitrary sign flips between neighboring surfaces.
    B = coeff(2);
    G = coeff(4);

    if B >= 0 && G >= 0
        transform = 'identity';
    elseif B < 0 && G < 0
        coeff([2, 4]) = -coeff([2, 4]);
        transform = 'theta-plus-pi';
    elseif B < 0
        coeff([2, 5]) = -coeff([2, 5]);
        transform = 'theta-reflection';
    else
        coeff([4, 5]) = -coeff([4, 5]);
        transform = 'pi-minus-theta';
    end
end

function theta = arc_length_angle(R, Z)
    Rclosed = [R(:); R(1)];
    Zclosed = [Z(:); Z(1)];
    ds = hypot(diff(Rclosed), diff(Zclosed));
    L = sum(ds);

    if ~isfinite(L) || L <= 0
        error('fit_CShape:DegenerateContour', ...
            'The input contour has zero or non-finite perimeter.');
    end

    s = [0; cumsum(ds(1:end-1))];
    theta = 2*pi*s/L;
end

function [R, Z] = prepare_contour(R, Z)
    R = R(:);
    Z = Z(:);

    if numel(R) ~= numel(Z) || any(~isfinite([R; Z]))
        error('fit_CShape:BadContour', ...
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
        error('fit_CShape:DuplicateContourPoint', ...
            'The contour contains coincident neighboring points.');
    end
end
