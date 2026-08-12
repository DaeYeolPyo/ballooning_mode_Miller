function [Gamma, theta_psi, diag] = evaluate_theta_psi( ...
        merluc, FFprime, pprime, qprime, opts)
    arguments
        merluc (1,1) struct
        FFprime (1,1) double
        pprime (1,1) double
        qprime (1,1) double
        opts.EnforcePeriodicity (1,1) logical = true
    end

    mu0 = 4*pi*1.e-7;
    theta_geo = merluc.theta_geo(:);
    theta_PEST = merluc.theta_PEST(:);

    F = merluc.F;
    q = merluc.q;
    R = merluc.R(:);
    Bp = merluc.Bp(:);

    % Gamma is affine in FFprime.  Separating the constant and FFprime
    % coefficients makes the Miller Eq. (21) periodicity constraint explicit.
    base = -qprime/q + ...
        (-2*merluc.invRc(:) - 2*merluc.sinu(:)./R + ...
         mu0*pprime*R./Bp)./(R.*Bp);
    FFprime_coefficient = 1/F^2 + 1./(R.^2.*Bp.^2);

    H = merluc.H(:);
    dl_dt = merluc.dl_dt(:);
    weight_geo = H.*dl_dt;

    periodic_denominator = trapz( ...
        theta_geo, weight_geo.*FFprime_coefficient);
    if abs(periodic_denominator) <= eps
        error('evaluate_theta_psi:SingularPeriodicityConstraint', ...
            'The Miller Eq. (21) FFprime coefficient is numerically zero.');
    end

    FFprime_periodic = -trapz(theta_geo, weight_geo.*base) ...
        / periodic_denominator;

    Gamma_input = base + FFprime.*FFprime_coefficient;
    if opts.EnforcePeriodicity
        FFprime_used = FFprime_periodic;
    else
        FFprime_used = FFprime;
    end

    Gamma = base + FFprime_used.*FFprime_coefficient;

    % Integrate H*Gamma*dl using the geometric parameter.  The equivalent
    % PEST-coordinate integral is integral Gamma*d(theta_PEST).
    theta_psi_input = cumtrapz( ...
        theta_geo, weight_geo.*Gamma_input);
    theta_psi = cumtrapz(theta_geo, weight_geo.*Gamma);
    theta_psi_PEST = cumtrapz(theta_PEST, Gamma);

    diag = struct();

    diag.theta = theta_PEST;
    diag.theta_geo = theta_geo;
    diag.theta_psi = theta_psi;
    diag.Gamma = Gamma;

    % Reference-surface geometry on the same closed [0,2*pi] grid.  The
    % normal basis is aligned with increasing Miller r/psi, while the tangent
    % basis follows increasing PEST theta.
    diag.R = R;
    diag.Z = merluc.Z(:);
    diag.Bp = Bp;
    diag.Bphi = merluc.Bphi(:);
    diag.B2 = merluc.B2(:);
    diag.B = merluc.B(:);
    diag.H = H;
    diag.dl_dt = dl_dt;
    diag.normal_R = merluc.sinu(:);
    diag.normal_Z = merluc.cosu(:);
    diag.tangent_R = merluc.Rt(:)./dl_dt;
    diag.tangent_Z = merluc.Zt(:)./dl_dt;
    diag.curvature_normal = merluc.invRc(:);
    diag.dBp_dl = periodic_derivative_closed(theta_geo, Bp)./dl_dt;

    diag.F = F;
    diag.q = q;
    diag.normal_sign = merluc.normal_sign;

    diag.FFprime_input = FFprime;
    diag.FFprime_periodic = FFprime_periodic;
    diag.FFprime_used = FFprime_used;
    diag.Gamma_input = Gamma_input;
    diag.theta_psi_input = theta_psi_input;
    diag.theta_psi_PEST = theta_psi_PEST;
    diag.closure_input = theta_psi_input(end) - theta_psi_input(1);
    diag.closure_used = theta_psi(end) - theta_psi(1);
    diag.closure_PEST = theta_psi_PEST(end) - theta_psi_PEST(1);
    diag.integration_max_abs_difference = max(abs( ...
        theta_psi - theta_psi_PEST));
    diag.input_closure_defect = abs(diag.closure_input) / max( ...
        trapz(theta_geo, abs(weight_geo.*Gamma_input)), eps);
    diag.periodicity_enforced = opts.EnforcePeriodicity;

    diag.FFprime = FFprime_used;
    diag.pprime = pprime;
    diag.qprime = qprime;
end

function df = periodic_derivative_closed(x, f)
    x = x(:);
    f = f(:);

    if numel(x) ~= numel(f) || numel(x) < 4
        error('evaluate_theta_psi:BadPeriodicDerivativeInput', ...
            'A closed periodic derivative requires at least four samples.');
    end

    period = x(end) - x(1);
    ncore = numel(x) - 1;
    h = period/ncore;
    fcore = f(1:ncore);
    dfcore = (circshift(fcore, -1) - circshift(fcore, 1))/(2*h);
    df = [dfcore; dfcore(1)];
end
