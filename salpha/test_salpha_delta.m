clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Equilibrium and scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_PT0.6');
%geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_circular');

target_psiN = 0.7;
delta_grid = [0.2, 0.4, 0.6];

% A 30 x 30 scan is intended to locate the marginal curves.  Increase the
% resolution, or refine only around lambda_max = 0, for final curves.
ntheta0 = 15;
nshat = 30;
shat_end = 5.0;
nalpha = 30;
alpha_end = 6.0;
theta_bnd = 5*pi;

n_surface_points = 1024;
n_balloon_points = 151;

%% Fit the reference surface once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);

surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=n_surface_points);
[param_reference, bnd_reference] = fit_Miller( ...
    eq, eqfunc, target_psiN, NTheta=n_surface_points);
merluc_reference = Miller_Mercier_Luc( ...
    surf, param_reference, bnd_reference);

FFprime = eqfunc.FFprime(target_psiN);
theta0 = linspace(-pi, pi, ntheta0 + 1);
theta0(end) = [];  % -pi and pi are the same ballooning phase
alpha = linspace(0, alpha_end, nalpha);
shat = linspace(0, shat_end, nshat);
theta_grid = linspace(-theta_bnd, theta_bnd, n_balloon_points).';

aN = eq.rmaxis;
BN = abs(eq.bcentr);
mu0 = 4*pi*1.e-7;

ndelta = numel(delta_grid);
lambda_max_delta = nan(ndelta, nshat, nalpha);
best_theta0_delta = nan(ndelta, nshat, nalpha);
lambda_theta0_delta = nan(ndelta, nshat, nalpha, ntheta0);
fail_message = strings(ndelta, nshat, nalpha, ntheta0);

volume_delta = nan(1, ndelta);
Vprime_delta = nan(1, ndelta);
dVdr_delta = nan(1, ndelta);
q_check_delta = nan(1, ndelta);
minimum_abs_jacobian = nan(1, ndelta);
param_delta = cell(1, ndelta);

nsolve = ndelta*nshat*nalpha*ntheta0;
fprintf('Reference Miller fit finished in %.2f s.\n', toc(setup_tic));
fprintf(['Reference delta = %.6g, kappa = %.6g, ', ...
    's_delta = %.6g.\n'], param_reference.delta, ...
    param_reference.kappa, param_reference.s_delta);
fprintf(['Starting delta s-alpha scan: %d delta x %d shat x ', ...
    '%d alpha x %d theta0 = %d solves.\n'], ...
    ndelta, nshat, nalpha, ntheta0, nsolve);

%% Delta -> shat -> alpha -> theta0 scan
scan_tic = tic;
total_rows = ndelta*nshat;

for dd = 1:ndelta
    delta_tic = tic;
    delta_value = delta_grid(dd);
    if ~isfinite(delta_value) || abs(delta_value) >= 1
        error('test_salpha_delta:BadDelta', ...
            'Every delta value must be finite and satisfy abs(delta) < 1.');
    end

    % Vary only delta.  The remaining fitted local-equilibrium parameters,
    % including s_delta = r*d(asin(delta))/dr, are held fixed.
    param_d = param_reference;
    param_d.delta = delta_value;
    bnd_d = Miller_boundary(param_d.R0, param_d.r, ...
        param_d.delta, param_d.kappa, NTheta=n_surface_points);
    merluc_d = Miller_Mercier_Luc(surf, param_d, bnd_d);
    param_delta{dd} = param_d;

    [volume_delta(dd), Vprime_delta(dd), dVdr_delta(dd)] = ...
        miller_volume_integrals(merluc_d);
    q_check_delta(dd) = merluc_d.q_check;
    minimum_abs_jacobian(dd) = min(abs(merluc_d.jac));

    if minimum_abs_jacobian(dd) <= ...
            1.e-10*max(abs(merluc_d.jac))
        error('test_salpha_delta:NearlySingularJacobian', ...
            'The Miller coordinate Jacobian is nearly singular at delta=%g.', ...
            delta_value);
    end

    fprintf(['\nDelta %d/%d: delta = %.4g, V = %.6g, ', ...
        'Vprime = %.6g, q_check/q = %.8f\n'], ...
        dd, ndelta, delta_value, volume_delta(dd), Vprime_delta(dd), ...
        q_check_delta(dd)/surf.q);

    for ss = 1:nshat
        row_tic = tic;
        qprime = shat(ss)*(surf.q*Vprime_delta(dd) ...
            /(2*volume_delta(dd)));

        for aa = 1:nalpha
            pprime = -alpha(aa)*(4*pi^2/(2*Vprime_delta(dd))) ...
                *sqrt(2*pi^2*eq.rmaxis/volume_delta(dd))/mu0;

            [~, ~, theta_diag] = evaluate_theta_psi(merluc_d, ...
                FFprime, pprime, qprime, EnforcePeriodicity=true);

            for tt = 1:ntheta0
                try
                    bal = evaluate_gcf(theta_diag, theta0(tt), ...
                        theta_bnd=theta_bnd, a_N=aN, B_N=BN);

                    % Resample the extended non-periodic coefficients onto a
                    % practical finite-element grid, as in test_salpha.m.
                    g = interp1(bal.theta, bal.g, theta_grid, 'pchip');
                    c = interp1(bal.theta, bal.c, theta_grid, 'pchip');
                    f = interp1(bal.theta, bal.f, theta_grid, 'pchip');

                    [Gmat, Cmat, Fmat, theta_Dof] = construct_matrix( ...
                        theta_grid, g, c, f);
                    [lambda, ~] = calculate_eigenvalue( ...
                        theta_Dof, Gmat, Cmat, Fmat, 'dirichlet');

                    lambda_theta0_delta(dd, ss, aa, tt) = lambda(1);
                catch ME
                    fail_message(dd, ss, aa, tt) = string(ME.message);
                end
            end

            values = reshape( ...
                lambda_theta0_delta(dd, ss, aa, :), 1, []);
            [lambda_max_delta(dd, ss, aa), imax] = ...
                max(values, [], 'omitnan');
            if isfinite(lambda_max_delta(dd, ss, aa))
                best_theta0_delta(dd, ss, aa) = theta0(imax);
            end
        end

        completed_rows = (dd - 1)*nshat + ss;
        elapsed = toc(scan_tic);
        eta = elapsed/completed_rows*(total_rows - completed_rows);
        fprintf(['  shat row %2d/%2d: shat = %.4g, ', ...
            'row %.1f s, total ETA %.1f min\n'], ...
            ss, nshat, shat(ss), toc(row_tic), eta/60);
    end

    fprintf('Delta %.4g scan finished in %.1f min.\n', ...
        delta_value, toc(delta_tic)/60);
end

elapsed = toc(scan_tic);
nfailed = nnz(strlength(fail_message) > 0);
fprintf('\nAll delta scans finished in %.1f min (%d failed solves).\n', ...
    elapsed/60, nfailed);

%% Package results
salpha_delta_result = struct();
salpha_delta_result.delta = delta_grid;
salpha_delta_result.shat = shat;
salpha_delta_result.alpha = alpha;
salpha_delta_result.theta0 = theta0;
salpha_delta_result.lambda_max = lambda_max_delta;
salpha_delta_result.lambda_theta0 = lambda_theta0_delta;
salpha_delta_result.best_theta0 = best_theta0_delta;
salpha_delta_result.fail_message = fail_message;
salpha_delta_result.volume = volume_delta;
salpha_delta_result.Vprime = Vprime_delta;
salpha_delta_result.dVdr = dVdr_delta;
salpha_delta_result.q_check = q_check_delta;
salpha_delta_result.minimum_abs_jacobian = minimum_abs_jacobian;
salpha_delta_result.reference_param = param_reference;
salpha_delta_result.param = param_delta;
salpha_delta_result.target_psiN = target_psiN;
salpha_delta_result.geqdsk_file = geqdsk_file;

%% Marginal-stability curves: max_theta0(lambda) = 0
figure('Color', 'w', 'Name', 'Delta s-alpha marginal stability curves');
ax = axes();
hold(ax, 'on');
colors = lines(ndelta);
curve_handles = gobjects(0);
curve_labels = strings(0);

for dd = 1:ndelta
    lambda_d = squeeze(lambda_max_delta(dd, :, :));
    has_crossing = min(lambda_d, [], 'all', 'omitnan') <= 0 ...
        && max(lambda_d, [], 'all', 'omitnan') >= 0;

    if has_crossing
        [~, h] = contour(ax, alpha, shat, lambda_d, [0, 0], ...
            'Color', colors(dd, :), 'LineWidth', 2.2);
        curve_handles(end + 1) = h; %#ok<SAGROW>
        curve_labels(end + 1) = ...
            sprintf('$\\delta = %.3g$', delta_grid(dd)); %#ok<SAGROW>
    else
        warning('test_salpha_delta:NoMarginalCrossing', ...
            ['No lambda_max=0 crossing was found for delta=%g. ', ...
             'Increase the shat or alpha range if necessary.'], ...
            delta_grid(dd));
    end
end

grid(ax, 'on');
box(ax, 'on');
ax.Layer = 'top';
ax.TickDir = 'out';
ax.FontSize = 11;
xlabel(ax, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax, '$\hat{s}$', 'Interpreter', 'latex');
if ~isempty(curve_handles)
    legend(ax, curve_handles, curve_labels, ...
        'Interpreter', 'latex', 'Location', 'best');
end

%% Local helper
function [V, Vprime, dVdr] = miller_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    % Enclosed toroidal volume and its local radial derivative:
    %   V      = pi*integral R^2*dZ
    %   dV/dr  = 2*pi*integral R*abs(J_{r,theta})*dtheta
    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVdr = 2*pi*trapz(theta, R.*abs(jac));

    dpsi_dr = abs(merluc.dpsi_dr);
    if ~isfinite(dpsi_dr) || dpsi_dr <= eps
        error('test_salpha_delta:BadFluxDerivative', ...
            'The Miller dpsi/dr value must be finite and nonzero.');
    end
    Vprime = dVdr/dpsi_dr;

    if any(~isfinite([V, Vprime, dVdr])) ...
            || any([V, Vprime, dVdr] <= 0)
        error('test_salpha_delta:BadVolumeIntegral', ...
            'The Miller V, Vprime, and dV/dr values must be positive.');
    end
end