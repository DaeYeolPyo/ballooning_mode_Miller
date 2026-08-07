function remapVal = validate_SFL_remap(surfaces, maps)
%VALIDATE_SFL_REMAP Check interpolation quality before/after SFL remap.

    names = ["PEST", "Hamada", "Boozer"];

    remapVal = struct();

    for n = 1:numel(names)
        name = names(n);
        map = maps.(name);

        Nrho = numel(surfaces);

        R_relerr      = nan(Nrho,1);
        Z_relerr      = nan(Nrho,1);
        lambda_relerr = nan(Nrho,1);

        R_maxerr      = nan(Nrho,1);
        Z_maxerr      = nan(Nrho,1);
        lambda_maxerr = nan(Nrho,1);

        theta_mono = false(Nrho,1);
        theta_uniform_err = theta_uniformity_error(map.theta);

        R_periodic_jump = nan(Nrho,1);
        Z_periodic_jump = nan(Nrho,1);
        lambda_periodic_jump = nan(Nrho,1);

        R_deriv_jump = nan(Nrho,1);
        Z_deriv_jump = nan(Nrho,1);
        lambda_deriv_jump = nan(Nrho,1);

        for k = 1:Nrho
            [thetaRaw, Rraw, Zraw, lambRaw] = get_raw_surface_data(surfaces(k), name);

            thetaRaw = thetaRaw(:);
            Rraw = Rraw(:);
            Zraw = Zraw(:);
            lambRaw = lambRaw(:);

            thetaRaw = thetaRaw - thetaRaw(1);
            thetaRaw = unwrap(thetaRaw);

            [thetaRaw, idx] = sort(thetaRaw);
            Rraw = Rraw(idx);
            Zraw = Zraw(idx);
            lambRaw = lambRaw(idx);

            keep = [true; diff(thetaRaw) > 1e-12];
            thetaRaw = thetaRaw(keep);
            Rraw = Rraw(keep);
            Zraw = Zraw(keep);
            lambRaw = lambRaw(keep);

            theta_mono(k) = all(diff(thetaRaw) > 0);

            Rmap = map.R(k,:).';
            Zmap = map.Z(k,:).';
            lmap = map.lambda(k,:).';

            % Compare remapped data back at the original raw theta locations.
            Rback = periodic_interp(map.theta(:), Rmap, thetaRaw);
            Zback = periodic_interp(map.theta(:), Zmap, thetaRaw);
            lback = periodic_interp(map.theta(:), lmap, thetaRaw);

            lambRaw = lambRaw - lambRaw(1);
            lback = lback - lback(1);

            Rerr = Rback - Rraw;
            Zerr = Zback - Zraw;
            lerr = lback - lambRaw;

            R_relerr(k) = rms(Rerr) / max(1, rms(Rraw));
            Z_relerr(k) = rms(Zerr) / max(1, rms(Zraw));
            lambda_relerr(k) = rms(lerr) / max(1, rms(lambRaw));

            R_maxerr(k) = max(abs(Rerr));
            Z_maxerr(k) = max(abs(Zerr));
            lambda_maxerr(k) = max(abs(lerr));

            per = periodic_endpoint_check(map.theta(:), Rmap, Zmap, lmap);

            R_periodic_jump(k) = per.R_jump;
            Z_periodic_jump(k) = per.Z_jump;
            lambda_periodic_jump(k) = per.lambda_jump;

            R_deriv_jump(k) = per.R_deriv_jump;
            Z_deriv_jump(k) = per.Z_deriv_jump;
            lambda_deriv_jump(k) = per.lambda_deriv_jump;
        end

        remapVal.(name) = struct();

        remapVal.(name).theta_monotone = theta_mono;
        remapVal.(name).theta_uniform_err = theta_uniform_err;

        remapVal.(name).R_relerr = R_relerr;
        remapVal.(name).Z_relerr = Z_relerr;
        remapVal.(name).lambda_relerr = lambda_relerr;

        remapVal.(name).R_maxerr = R_maxerr;
        remapVal.(name).Z_maxerr = Z_maxerr;
        remapVal.(name).lambda_maxerr = lambda_maxerr;

        remapVal.(name).R_periodic_jump = R_periodic_jump;
        remapVal.(name).Z_periodic_jump = Z_periodic_jump;
        remapVal.(name).lambda_periodic_jump = lambda_periodic_jump;

        remapVal.(name).R_deriv_jump = R_deriv_jump;
        remapVal.(name).Z_deriv_jump = Z_deriv_jump;
        remapVal.(name).lambda_deriv_jump = lambda_deriv_jump;
    end

    print_remap_validation_summary(remapVal);
end

function yq = periodic_interp(theta, y, thetaq)
    theta = theta(:);
    y = y(:);
    thetaq = thetaq(:);

    theta = theta - theta(1);
    theta = unwrap(theta);

    theta_ext = [theta; 2*pi];
    y_ext = [y; y(1)];

    thetaq = mod(thetaq, 2*pi);

    yq = interp1(theta_ext, y_ext, thetaq, 'pchip');
end

function err = theta_uniformity_error(theta)
    theta = theta(:);
    dtheta = diff([theta; theta(1) + 2*pi]);
    err = max(abs(dtheta - mean(dtheta))) / max(eps, abs(mean(dtheta)));
end

function per = periodic_endpoint_check(theta, R, Z, lambda)
    theta = theta(:);
    R = R(:);
    Z = Z(:);
    lambda = lambda(:);

    dtheta = theta(2) - theta(1);

    % Extrapolated endpoint at 2*pi from last grid point.
    R_2pi = periodic_interp(theta, R, 2*pi);
    Z_2pi = periodic_interp(theta, Z, 2*pi);
    l_2pi = periodic_interp(theta, lambda, 2*pi);

    per.R_jump = abs(R_2pi - R(1)) / max(1, max(abs(R)));
    per.Z_jump = abs(Z_2pi - Z(1)) / max(1, max(abs(Z)));
    per.lambda_jump = abs(l_2pi - lambda(1)) / max(1, max(abs(lambda)));

    dR = periodic_derivative_uniform(R, dtheta);
    dZ = periodic_derivative_uniform(Z, dtheta);
    dl = periodic_derivative_uniform(lambda, dtheta);

    per.R_deriv_jump = abs(dR(1) - dR(end)) / max(1, max(abs(dR)));
    per.Z_deriv_jump = abs(dZ(1) - dZ(end)) / max(1, max(abs(dZ)));
    per.lambda_deriv_jump = abs(dl(1) - dl(end)) / max(1, max(abs(dl)));
end

function dydt = periodic_derivative_uniform(y, dtheta)
    y = y(:);
    dydt = (circshift(y, -1) - circshift(y, 1)) ./ (2*dtheta);
end

function print_remap_validation_summary(remapVal)
    names = string(fieldnames(remapVal));

    fprintf('\n===== SFL remap validation =====\n');
    fprintf('%10s %12s %12s %12s %12s %12s %12s\n', ...
        'coord', 'R rms', 'Z rms', 'lam rms', ...
        'R per', 'Z per', 'lam per');

    for n = 1:numel(names)
        name = names(n);
        v = remapVal.(name);

        fprintf('%10s %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n', ...
            name, ...
            median(v.R_relerr, 'omitnan'), ...
            median(v.Z_relerr, 'omitnan'), ...
            median(v.lambda_relerr, 'omitnan'), ...
            median(v.R_periodic_jump, 'omitnan'), ...
            median(v.Z_periodic_jump, 'omitnan'), ...
            median(v.lambda_periodic_jump, 'omitnan'));
    end
end