clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% NT reference equilibrium and elongation scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_NT0.6');
target_psiN = 0.8;
miller_fit_dpsi = 8.e-2;
kappa_grid = [1.3, 1.5, 1.7];

ntheta0 = 16;
nshat = 20;
shat_end = 7.0;
nalpha = 20;
alpha_end = 12.0;
theta_bnd = 5*pi;
n_surface_points = int32(1024);
n_balloon_points = 151;

%% Fit the NT Miller reference surface once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=n_surface_points);

[param_reference, bnd_reference] = fit_Miller( ...
    eq, eqfunc, target_psiN, NTheta=n_surface_points, ...
    dpsi=miller_fit_dpsi);
merluc_reference = Miller_Mercier_Luc( ...
    surf, param_reference, bnd_reference);

eq_on_reference = eqfunc.eval( ...
    merluc_reference.R, merluc_reference.Z);
Bp_reference_relative_rms = relative_rms( ...
    merluc_reference.Bp, eq_on_reference.Bp);

FFprime = eqfunc.FFprime(target_psiN);
qprime_equilibrium = eqfunc.qprime(target_psiN);
pprime_equilibrium = eqfunc.pprime(target_psiN);

theta0 = 2*pi*((0:ntheta0-1) - floor(ntheta0/2))/ntheta0;
alpha = linspace(0, alpha_end, nalpha);
shat = linspace(0, shat_end, nshat);
theta_grid = linspace(-theta_bnd, theta_bnd, n_balloon_points).';
aN = eq.rmaxis;
BN = abs(eq.bcentr);
mu0 = 4*pi*1.e-7;

nkappa = numel(kappa_grid);
lambda_max_kappa = nan(nkappa, nshat, nalpha);
best_theta0_kappa = nan(nkappa, nshat, nalpha);
lambda_theta0_kappa = nan(nkappa, nshat, nalpha, ntheta0);
fail_message = strings(nkappa, nshat, nalpha, ntheta0);

volume_kappa = nan(1, nkappa);
Vprime_kappa = nan(1, nkappa);
dVdr_kappa = nan(1, nkappa);
dpsi_dr_kappa = nan(1, nkappa);
q_check_kappa = nan(1, nkappa);
minimum_abs_jacobian = nan(1, nkappa);
qprime_scan = nan(nkappa, nshat);
pprime_scan = nan(nkappa, nalpha);
shat_equilibrium = nan(1, nkappa);
alpha_equilibrium = nan(1, nkappa);
param_kappa = cell(1, nkappa);
merluc_kappa = cell(1, nkappa);

fprintf('NT Miller kappa scan reference\n');
fprintf('  GEQDSK                      : %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  reference delta / kappa     : %+.8f / %.8f\n', ...
    param_reference.delta, param_reference.kappa);
fprintf('  reference s_delta / s_kappa : %+.8f / %+.8f\n', ...
    param_reference.s_delta, param_reference.s_kappa);
fprintf('  reference Bp relative RMS   : %.6e\n', ...
    Bp_reference_relative_rms);
fprintf('  radial-fit half width       : %.6f in psiN\n', ...
    param_reference.radial_fit.half_width);
fprintf('  scan kappa values           :');
fprintf(' %.1f', kappa_grid);
fprintf('\n');

assert(param_reference.delta < 0, ...
    'The fitted reference surface is not negative triangularity.');
assert(Bp_reference_relative_rms < 2.e-2, ...
    'The reference NT Miller Bp fit exceeds 2%% RMS error.');

%% Build all modified NT local equilibria
for kk = 1:nkappa
    kappa = kappa_grid(kk);
    param_k = param_reference;
    param_k.kappa = kappa;

    % Miller_Mercier_Luc uses s_kappa=(r/kappa)*dkappa/dr, so retaining
    % the fitted s_kappa automatically updates dkappa/dr consistently.
    bnd_k = Miller_boundary( ...
        param_k.R0, param_k.r, param_k.delta, param_k.kappa, ...
        Ntheta=n_surface_points);
    merluc_k = Miller_Mercier_Luc(surf, param_k, bnd_k);

    [volume_kappa(kk), Vprime_kappa(kk), dVdr_kappa(kk)] = ...
        miller_volume_integrals(merluc_k);
    dpsi_dr_kappa(kk) = merluc_k.dpsi_dr;
    q_check_kappa(kk) = merluc_k.q_check;
    minimum_abs_jacobian(kk) = min(abs(merluc_k.jac));

    normal_projection_r = ...
        merluc_k.normal_sign*merluc_k.normal_projection_raw;
    normal_projection_psi = ...
        sign(merluc_k.dpsi_dr)*normal_projection_r;
    if any(normal_projection_psi <= 0)
        error('test_salpha_NT_kappa:BadNormal', ...
            'The Mercier normal is misaligned at kappa=%g.', kappa);
    end

    jacobian_scale = max(abs(merluc_k.jac));
    if minimum_abs_jacobian(kk) <= 1.e-10*jacobian_scale
        error('test_salpha_NT_kappa:SingularJacobian', ...
            'The Miller coordinate is nearly singular at kappa=%g.', ...
            kappa);
    end
    if abs(q_check_kappa(kk) - surf.q) ...
            > 1.e-10*max(1, abs(surf.q))
        error('test_salpha_NT_kappa:BadQ', ...
            'The modified NT surface does not reproduce q at kappa=%g.', ...
            kappa);
    end

    qprime_scan(kk, :) = shat.*( ...
        surf.q*Vprime_kappa(kk)/(2*volume_kappa(kk)));
    pprime_scan(kk, :) = -alpha.*(4*pi^2/(2*Vprime_kappa(kk))) ...
        .*sqrt(2*pi^2*eq.rmaxis/volume_kappa(kk))/mu0;
    shat_equilibrium(kk) = 2*volume_kappa(kk)*qprime_equilibrium ...
        /(surf.q*Vprime_kappa(kk));
    alpha_equilibrium(kk) = -mu0*pprime_equilibrium ...
        *(2*Vprime_kappa(kk))/(4*pi^2) ...
        *sqrt(volume_kappa(kk)/(2*pi^2*eq.rmaxis));

    param_kappa{kk} = param_k;
    merluc_kappa{kk} = merluc_k;

    fprintf(['  kappa=%.1f: V=%.6e, Vprime=%+.6e, ', ...
        'dpsi/dr=%+.6e, min|J|=%.3e\n'], ...
        kappa, volume_kappa(kk), Vprime_kappa(kk), ...
        dpsi_dr_kappa(kk), minimum_abs_jacobian(kk));
end

%% kappa -> shat -> alpha -> theta0 scan
nsolve = nkappa*nshat*nalpha*ntheta0;
fprintf('Reference setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting NT kappa scan: %d kappa x %d shat x %d alpha x ', ...
    '%d theta0 = %d solves.\n'], ...
    nkappa, nshat, nalpha, ntheta0, nsolve);

scan_tic = tic;
total_rows = nkappa*nshat;

for kk = 1:nkappa
    kappa_tic = tic;
    merluc_k = merluc_kappa{kk};
    fprintf('\nKappa %d/%d: %.1f\n', kk, nkappa, kappa_grid(kk));

    for ss = 1:nshat
        row_tic = tic;
        qprime = qprime_scan(kk, ss);

        for aa = 1:nalpha
            pprime = pprime_scan(kk, aa);
            [~, ~, theta_diag] = evaluate_theta_psi(merluc_k, ...
                FFprime, pprime, qprime, EnforcePeriodicity=true);

            for tt = 1:ntheta0
                try
                    bal = evaluate_gcf(theta_diag, theta0(tt), ...
                        theta_bnd=theta_bnd, a_N=aN, B_N=BN);
                    g = interp1( ...
                        bal.theta, bal.g, theta_grid, 'pchip');
                    c = interp1( ...
                        bal.theta, bal.c, theta_grid, 'pchip');
                    f = interp1( ...
                        bal.theta, bal.f, theta_grid, 'pchip');

                    [Gmat, Cmat, Fmat, theta_Dof] = ...
                        construct_matrix(theta_grid, g, c, f);
                    [lambda, ~] = calculate_eigenvalue( ...
                        theta_Dof, Gmat, Cmat, Fmat, 'dirichlet');
                    lambda_theta0_kappa(kk, ss, aa, tt) = lambda(1);
                catch ME
                    fail_message(kk, ss, aa, tt) = string(ME.message);
                end
            end

            values = reshape( ...
                lambda_theta0_kappa(kk, ss, aa, :), 1, []);
            [lambda_max_kappa(kk, ss, aa), imax] = ...
                max(values, [], 'omitnan');
            if isfinite(lambda_max_kappa(kk, ss, aa))
                best_theta0_kappa(kk, ss, aa) = theta0(imax);
            end
        end

        completed_rows = (kk - 1)*nshat + ss;
        elapsed = toc(scan_tic);
        eta = elapsed/completed_rows*(total_rows - completed_rows);
        fprintf(['  shat row %2d/%2d: shat=%.4g, row %.1f s, ', ...
            'total ETA %.1f min\n'], ss, nshat, shat(ss), ...
            toc(row_tic), eta/60);
    end

    fprintf('Kappa=%.1f finished in %.1f min.\n', ...
        kappa_grid(kk), toc(kappa_tic)/60);
end

elapsed = toc(scan_tic);
nfailed = nnz(strlength(fail_message) > 0);
fprintf(['\nAll NT kappa scans finished in %.1f min ', ...
    '(%d failed solves).\n'], elapsed/60, nfailed);

%% Package and validate results
salpha_NT_kappa_result = struct();
salpha_NT_kappa_result.model = ...
    'negative-triangularity-kappa-scan-miller-mercier-luc';
salpha_NT_kappa_result.kappa = kappa_grid;
salpha_NT_kappa_result.s_kappa = param_reference.s_kappa;
salpha_NT_kappa_result.shat = shat;
salpha_NT_kappa_result.alpha = alpha;
salpha_NT_kappa_result.theta0 = theta0;
salpha_NT_kappa_result.lambda_max = lambda_max_kappa;
salpha_NT_kappa_result.lambda_theta0 = lambda_theta0_kappa;
salpha_NT_kappa_result.best_theta0 = best_theta0_kappa;
salpha_NT_kappa_result.fail_message = fail_message;
salpha_NT_kappa_result.qprime = qprime_scan;
salpha_NT_kappa_result.pprime = pprime_scan;
salpha_NT_kappa_result.volume = volume_kappa;
salpha_NT_kappa_result.Vprime = Vprime_kappa;
salpha_NT_kappa_result.dVdr = dVdr_kappa;
salpha_NT_kappa_result.dpsi_dr = dpsi_dr_kappa;
salpha_NT_kappa_result.q_check = q_check_kappa;
salpha_NT_kappa_result.minimum_abs_jacobian = ...
    minimum_abs_jacobian;
salpha_NT_kappa_result.shat_equilibrium = shat_equilibrium;
salpha_NT_kappa_result.alpha_equilibrium = alpha_equilibrium;
salpha_NT_kappa_result.reference_param = param_reference;
salpha_NT_kappa_result.param = param_kappa;
salpha_NT_kappa_result.target_psiN = target_psiN;
salpha_NT_kappa_result.geqdsk_file = geqdsk_file;
salpha_NT_kappa_result.Bp_reference_relative_rms = ...
    Bp_reference_relative_rms;
salpha_NT_kappa_result.scan_elapsed_seconds = elapsed;
salpha_NT_kappa_result.nfailed = nfailed;

has_marginal_crossing = false(1, nkappa);
lambda_range = nan(nkappa, 2);
for kk = 1:nkappa
    lambda_k = squeeze(lambda_max_kappa(kk, :, :));
    finite_k = lambda_k(isfinite(lambda_k));
    if ~isempty(finite_k)
        lambda_range(kk, :) = [min(finite_k), max(finite_k)];
        has_marginal_crossing(kk) = lambda_range(kk, 1) <= 0 ...
            && lambda_range(kk, 2) >= 0;
    end
end
salpha_NT_kappa_result.lambda_range = lambda_range;
salpha_NT_kappa_result.has_marginal_crossing = ...
    has_marginal_crossing;

assert(nfailed == 0, ...
    'The NT kappa scan contains failed eigenvalue solves.');
assert(all(isfinite(lambda_max_kappa), 'all'), ...
    'The NT kappa scan contains non-finite eigenvalues.');
assert(all(lambda_max_kappa(:, :, 1) < 0, 'all'), ...
    'The alpha=0 column should be stable for every kappa.');

fprintf('\nNT kappa scan summary\n');
for kk = 1:nkappa
    fprintf(['  kappa=%.1f: lambda=[%+.4e,%+.4e], crossing=%d, ', ...
        '(shat,alpha)_eq=(%.4f,%.4f)\n'], ...
        kappa_grid(kk), lambda_range(kk, 1), lambda_range(kk, 2), ...
        has_marginal_crossing(kk), ...
        shat_equilibrium(kk), alpha_equilibrium(kk));
end

%% Individual stability maps with a shared color scale
finite_lambda = lambda_max_kappa(isfinite(lambda_max_kappa));
color_limits = [min(finite_lambda), max(finite_lambda)];
map_figure = figure('Color', 'w', 'Name', ...
    'NT kappa scan stability maps');
map_layout = tiledlayout(map_figure, 1, nkappa, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
map_axes = gobjects(1, nkappa);

for kk = 1:nkappa
    map_axes(kk) = nexttile(map_layout);
    lambda_k = squeeze(lambda_max_kappa(kk, :, :));
    contourf(map_axes(kk), alpha, shat, lambda_k, 30, ...
        'LineStyle', 'none');
    hold(map_axes(kk), 'on');
    if has_marginal_crossing(kk)
        contour(map_axes(kk), alpha, shat, lambda_k, [0, 0], ...
            'k', 'LineWidth', 1.8);
    end
    plot(map_axes(kk), alpha_equilibrium(kk), ...
        shat_equilibrium(kk), 'p', ...
        'MarkerSize', 9, 'MarkerEdgeColor', 'w', ...
        'MarkerFaceColor', [0.85, 0.1, 0.75], 'LineWidth', 1.0);
    clim(map_axes(kk), color_limits);
    grid(map_axes(kk), 'on');
    box(map_axes(kk), 'on');
    map_axes(kk).Layer = 'top';
    xlabel(map_axes(kk), '$\alpha$', 'Interpreter', 'latex');
    ylabel(map_axes(kk), '$\hat{s}$', 'Interpreter', 'latex');
    title(map_axes(kk), sprintf('$\\kappa=%.1f$', kappa_grid(kk)), ...
        'Interpreter', 'latex');
end

cb_map = colorbar(map_axes(end));
cb_map.Label.Interpreter = 'latex';
cb_map.Label.String = '$\lambda_{\max}$';
title(map_layout, sprintf(['NT Miller elongation scan, ', ...
    '$\\psi_N=%.2f$'], target_psiN), 'Interpreter', 'latex');

%% Overlay marginal-stability contours
marginal_figure = figure('Color', 'w', 'Name', ...
    'NT kappa marginal curves');
ax_marginal = axes(marginal_figure);
hold(ax_marginal, 'on');
colors = lines(nkappa);
curve_handles = gobjects(0);
curve_labels = strings(0);

for kk = 1:nkappa
    lambda_k = squeeze(lambda_max_kappa(kk, :, :));
    if has_marginal_crossing(kk)
        [~, curve] = contour(ax_marginal, alpha, shat, ...
            lambda_k, [0, 0], 'Color', colors(kk, :), ...
            'LineWidth', 2.2);
        curve_handles(end + 1) = curve; %#ok<SAGROW>
        curve_labels(end + 1) = ...
            sprintf('$\\kappa=%.1f$', kappa_grid(kk)); %#ok<SAGROW>
    else
        warning('test_salpha_NT_kappa:NoCrossing', ...
            ['No marginal crossing was found for kappa=%g. ', ...
             'Increase the scan range.'], kappa_grid(kk));
    end
end

grid(ax_marginal, 'on');
box(ax_marginal, 'on');
ax_marginal.Layer = 'top';
ax_marginal.TickDir = 'out';
xlabel(ax_marginal, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax_marginal, '$\hat{s}$', 'Interpreter', 'latex');
title(ax_marginal, sprintf(['NT Miller elongation scan, ', ...
    '$\\psi_N=%.2f$, fixed $s_\\kappa=%.3f$'], ...
    target_psiN, param_reference.s_kappa), 'Interpreter', 'latex');
if ~isempty(curve_handles)
    legend(ax_marginal, curve_handles, curve_labels, ...
        'Interpreter', 'latex', 'Location', 'best');
end

salpha_NT_kappa_result.map_figure = map_figure;
salpha_NT_kappa_result.marginal_figure = marginal_figure;
salpha_NT_kappa_result.status = 'PASS';

fprintf('  validation status           : PASS\n');

function [V, Vprime, dVdr] = miller_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVdr = 2*pi*trapz(theta, R.*abs(jac));

    if ~isfinite(merluc.dpsi_dr) || abs(merluc.dpsi_dr) <= eps
        error('test_salpha_NT_kappa:BadFluxDerivative', ...
            'The signed Miller dpsi/dr must be finite and nonzero.');
    end
    Vprime = dVdr/merluc.dpsi_dr;

    if any(~isfinite([V, Vprime, dVdr])) ...
            || V <= 0 || dVdr <= 0 || Vprime == 0
        error('test_salpha_NT_kappa:BadVolumeIntegral', ...
            'A modified V, Vprime, or dV/dr value is invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
