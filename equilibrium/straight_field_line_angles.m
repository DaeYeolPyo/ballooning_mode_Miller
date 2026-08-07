function ang = straight_field_line_angles(eqfunc, surf, ints, opts)
%COMPUTE_SFL_ANGLES Compute PEST, Hamada, and Boozer angles on one surface.
%
% Inputs:
%   I    : build_eq_interpolants output
%   fs   : extract_flux_surface output
%   ints : compute_q_V_Vprime output
%
% Output:
%   ang.theta_geom
%   ang.theta_PEST
%   ang.theta_Hamada
%   ang.theta_Boozer
%   ang.lambda_Hamada
%   ang.lambda_Boozer

    arguments
        eqfunc struct
        surf struct
        ints struct
        opts.UseQ string = "integral"   % "integral" or "geqdsk"
    end

    R  = surf.R(:);
    Z  = surf.Z(:);
    dl = surf.dl(:);

    val = eqfunc.eval(R, Z);

    gradPsi = val.gradPsi;
    F       = val.F;
    B2      = val.B2;

    if any(~isfinite(gradPsi)) || any(gradPsi <= 0)
        error('compute_sfl_angles:BadGradPsi', ...
              'Invalid |grad psi| found on flux surface.');
    end

    wPEST   = F ./ (R .* gradPsi);
    wHamada = R ./ gradPsi;
    wBoozer = R .* B2 ./ gradPsi;

    theta_PEST   = normalized_cumulative_angle(wPEST,   dl);
    theta_Hamada = normalized_cumulative_angle(wHamada, dl);
    theta_Boozer = normalized_cumulative_angle(wBoozer, dl);

    theta_geom = geometric_poloidal_angle(R, Z, surf.Raxis, surf.Zaxis);

    switch opts.UseQ
        case "integral"
            q = ints.q_int;
        case "geqdsk"
            q = ints.q_geqdsk;
        otherwise
            error('compute_sfl_angles:BadUseQ', ...
                  'opts.UseQ must be "integral" or "geqdsk".');
    end

    lambda_Hamada = q * unwrap_periodic_difference(theta_Hamada - theta_PEST);
    lambda_Boozer = q * unwrap_periodic_difference(theta_Boozer - theta_PEST);

    % Remove arbitrary flux-function offset so lambda starts at zero.
    lambda_Hamada = lambda_Hamada - lambda_Hamada(1);
    lambda_Boozer = lambda_Boozer - lambda_Boozer(1);

    ang = struct();

    ang.psi  = surf.psi;
    ang.psiN = surf.psiN;

    ang.theta_geom   = theta_geom;
    ang.theta_PEST   = theta_PEST;
    ang.theta_Hamada = theta_Hamada;
    ang.theta_Boozer = theta_Boozer;

    ang.lambda_Hamada = lambda_Hamada;
    ang.lambda_Boozer = lambda_Boozer;

    ang.q_used = q;

    ang.weights = struct();
    ang.weights.PEST   = wPEST;
    ang.weights.Hamada = wHamada;
    ang.weights.Boozer = wBoozer;

    ang.diagnostics = struct();
    ang.diagnostics.PEST_monotone   = all(diff(theta_PEST)   > 0);
    ang.diagnostics.Hamada_monotone = all(diff(theta_Hamada) > 0);
    ang.diagnostics.Boozer_monotone = all(diff(theta_Boozer) > 0);
end

function theta = normalized_cumulative_angle(w, dl)
%NORMALIZED_CUMULATIVE_ANGLE Build theta in [0, 2*pi) from contour weights.
%
% w is defined at contour points.
% dl(i) is the segment length from point i to point i+1 periodically.

    w  = w(:);
    dl = dl(:);

    wp = circshift(w, -1);

    % Segment-centered trapezoidal weight.
    dA = 0.5 * (w + wp) .* dl;

    total = sum(dA);

    if ~isfinite(total) || abs(total) < eps
        error('compute_sfl_angles:ZeroWeightIntegral', ...
              'Weight integral is zero or invalid.');
    end

    cum = [0; cumsum(dA(1:end-1))];

    theta = 2*pi * cum / total;

    % If both cumulative and total are negative this is already increasing.
    % Normalize small roundoff.
    theta = mod(theta, 2*pi);
    theta(1) = 0;

    theta = unwrap(theta);
    theta = theta - theta(1);
end

function theta = geometric_poloidal_angle(R, Z, Raxis, Zaxis)
%GEOMETRIC_POLOIDAL_ANGLE Geometric angle around magnetic axis.

    theta = atan2(Z - Zaxis, R - Raxis);
    theta = unwrap(theta);

    theta = theta - theta(1);

    if theta(end) < 0
        theta = theta + 2*pi;
    end

    theta = mod(theta, 2*pi);
    theta(1) = 0;
    theta = unwrap(theta);
end

function x = unwrap_periodic_difference(x)
%UNWRAP_PERIODIC_DIFFERENCE Keep angle differences smooth.

    x = unwrap(x(:));
    x = x - x(1);
end