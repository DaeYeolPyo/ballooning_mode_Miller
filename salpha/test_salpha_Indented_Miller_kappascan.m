clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Reference equilibrium and elongation scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_I0.8');
target_psiN = 0.7;
kappa_grid = [1.3, 1.5, 1.7, 1.9];

ntheta0 = 16;
nshat = 30;
shat_end = 7.0;
nalpha = 30;
alpha_end = 12.0;
theta_bnd = 5*pi;
n_surface_points = int32(1024);
n_fit_points = int32(256);
n_balloon_points = 151;

%% Fit the reference indented-Miller equilibrium once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=n_surface_points);

[param_reference, bnd_reference] = fit_Indented_Miller( ...
    eq, eqfunc, target_psiN, ...
    NTheta=n_fit_points, NRadialFit=int32(9), dpsi=3.e-2, ...
    MaxIterations=int32(250), NPhaseStarts=int32(32), ...
    NXStarts=int32(33));
merluc_reference = Indented_Miller_Mercier_Luc( ...
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
dVda_kappa = nan(1, nkappa);
dpsi_da_kappa = nan(1, nkappa);
q_check_kappa = nan(1, nkappa);
minimum_abs_jacobian = nan(1, nkappa);
qprime_scan = nan(nkappa, nshat);
pprime_scan = nan(nkappa, nalpha);
shat_equilibrium = nan(1, nkappa);
alpha_equilibrium = nan(1, nkappa);
param_kappa = cell(1, nkappa);
merluc_kappa = cell(1, nkappa);

fprintf('Indented-Miller kappa scan reference\n');
fprintf('  GEQDSK                      : %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  reference kappa             : %.8f\n', ...
    param_reference.kappa);
fprintf('  reference s_kappa           : %+.8f\n', ...
    param_reference.s_kappa);
fprintf('  reference delta / I0        : %+.8f / %.8f\n', ...
    param_reference.delta, param_reference.I0);
fprintf('  reference Bp relative RMS   : %.6e\n', ...
    Bp_reference_relative_rms);
fprintf('  scan kappa values           :');
fprintf(' %.1f', kappa_grid);
fprintf('\n');

assert(Bp_reference_relative_rms < 2.e-2, ...
    'The reference indented-Miller Bp fit exceeds 2%% RMS error.');

%% Build all modified local equilibria
for kk = 1:nkappa
    kappa = kappa_grid(kk);
    param_k = set_elongation( ...
        param_reference, kappa, param_reference.s_kappa);
    bnd_k = Indented_Miller_boundary( ...
        param_k.R0, param_k.a, param_k.delta, ...
        param_k.kappa, param_k.I0, NTheta=n_fit_points);
    merluc_k = Indented_Miller_Mercier_Luc(surf, param_k, bnd_k);

    [volume_kappa(kk), Vprime_kappa(kk), dVda_kappa(kk)] = ...
        indented_volume_integrals(merluc_k);
    dpsi_da_kappa(kk) = merluc_k.dpsi_da;
    q_check_kappa(kk) = merluc_k.q_check;
    minimum_abs_jacobian(kk) = min(abs(merluc_k.jac));

    normal_projection_a = ...
        merluc_k.normal_sign*merluc_k.normal_projection_raw;
    normal_projection_psi = ...
        sign(merluc_k.dpsi_da)*normal_projection_a;
    if any(normal_projection_psi <= 0)
        error('test_salpha_Indented_Miller_kappascan:BadNormal', ...
            'The Mercier normal is misaligned at kappa=%g.', kappa);
    end

    jacobian_scale = max(abs(merluc_k.jac));
    if minimum_abs_jacobian(kk) <= 1.e-10*jacobian_scale
        error('test_salpha_Indented_Miller_kappascan:SingularJacobian', ...
            'The fitted coordinate is nearly singular at kappa=%g.', ...
            kappa);
    end
    if abs(q_check_kappa(kk) - surf.q) ...
            > 1.e-10*max(1, abs(surf.q))
        error('test_salpha_Indented_Miller_kappascan:BadQ', ...
            'The modified surface does not reproduce q at kappa=%g.', ...
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
        'dpsi/da=%+.6e, min|J|=%.3e\n'], ...
        kappa, volume_kappa(kk), Vprime_kappa(kk), ...
        dpsi_da_kappa(kk), minimum_abs_jacobian(kk));
end

%% kappa -> shat -> alpha -> theta0 scan
nsolve = nkappa*nshat*nalpha*ntheta0;
fprintf('Reference setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting kappa scan: %d kappa x %d shat x %d alpha x ', ...
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
fprintf(['\nAll kappa scans finished in %.1f min ', ...
    '(%d failed solves).\n'], elapsed/60, nfailed);

%% Package and validate results
salpha_Indented_Miller_kappascan_result = struct();
salpha_Indented_Miller_kappascan_result.model = ...
    'indented-miller-kappa-scan-mercier-luc';
salpha_Indented_Miller_kappascan_result.kappa = kappa_grid;
salpha_Indented_Miller_kappascan_result.s_kappa = ...
    param_reference.s_kappa;
salpha_Indented_Miller_kappascan_result.shat = shat;
salpha_Indented_Miller_kappascan_result.alpha = alpha;
salpha_Indented_Miller_kappascan_result.theta0 = theta0;
salpha_Indented_Miller_kappascan_result.lambda_max = ...
    lambda_max_kappa;
salpha_Indented_Miller_kappascan_result.lambda_theta0 = ...
    lambda_theta0_kappa;
salpha_Indented_Miller_kappascan_result.best_theta0 = ...
    best_theta0_kappa;
salpha_Indented_Miller_kappascan_result.fail_message = fail_message;
salpha_Indented_Miller_kappascan_result.qprime = qprime_scan;
salpha_Indented_Miller_kappascan_result.pprime = pprime_scan;
salpha_Indented_Miller_kappascan_result.volume = volume_kappa;
salpha_Indented_Miller_kappascan_result.Vprime = Vprime_kappa;
salpha_Indented_Miller_kappascan_result.dVda = dVda_kappa;
salpha_Indented_Miller_kappascan_result.dpsi_da = dpsi_da_kappa;
salpha_Indented_Miller_kappascan_result.q_check = q_check_kappa;
salpha_Indented_Miller_kappascan_result.minimum_abs_jacobian = ...
    minimum_abs_jacobian;
salpha_Indented_Miller_kappascan_result.shat_equilibrium = ...
    shat_equilibrium;
salpha_Indented_Miller_kappascan_result.alpha_equilibrium = ...
    alpha_equilibrium;
salpha_Indented_Miller_kappascan_result.reference_param = ...
    param_reference;
salpha_Indented_Miller_kappascan_result.param = param_kappa;
salpha_Indented_Miller_kappascan_result.target_psiN = target_psiN;
salpha_Indented_Miller_kappascan_result.geqdsk_file = geqdsk_file;
salpha_Indented_Miller_kappascan_result.Bp_reference_relative_rms = ...
    Bp_reference_relative_rms;
salpha_Indented_Miller_kappascan_result.scan_elapsed_seconds = elapsed;
salpha_Indented_Miller_kappascan_result.nfailed = nfailed;

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
salpha_Indented_Miller_kappascan_result.lambda_range = lambda_range;
salpha_Indented_Miller_kappascan_result.has_marginal_crossing = ...
    has_marginal_crossing;

assert(nfailed == 0, ...
    'The kappa scan contains failed eigenvalue solves.');
assert(all(isfinite(lambda_max_kappa), 'all'), ...
    'The kappa scan contains non-finite eigenvalues.');
assert(all(lambda_max_kappa(:, :, 1) < 0, 'all'), ...
    'The alpha=0 column should be stable for every kappa.');

fprintf('\nKappa scan summary\n');
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
    'Indented-Miller kappa scan stability maps');
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

%% Overlay marginal-stability contours
marginal_figure = figure('Color', 'w', 'Name', ...
    'Indented-Miller kappa marginal curves');
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
        warning('test_salpha_Indented_Miller_kappascan:NoCrossing', ...
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

if ~isempty(curve_handles)
    legend(ax_marginal, curve_handles, curve_labels, ...
        'Interpreter', 'latex', 'Location', 'best');
end

salpha_Indented_Miller_kappascan_result.map_figure = map_figure;
salpha_Indented_Miller_kappascan_result.marginal_figure = ...
    marginal_figure;
salpha_Indented_Miller_kappascan_result.status = 'PASS';

fprintf('  validation status           : PASS\n');

function param = set_elongation(param, kappa, s_kappa)
    param.kappa = kappa;
    param.b = param.a*kappa;
    param.coeff(4) = param.b;

    % Preserve s_kappa=(a/kappa)*dkappa/da while changing kappa.
    param.s_kappa = s_kappa;
    param.dkappa_da = s_kappa*kappa/param.a;
    param.db_da = kappa*(1 + s_kappa);
    param.dcoeff_da(4) = param.db_da;
    param.coeff_table = [param.coeff, param.dcoeff_da];
end

function [V, Vprime, dVda] = indented_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVda = 2*pi*trapz(theta, R.*abs(jac));

    if ~isfinite(merluc.dpsi_da) || abs(merluc.dpsi_da) <= eps
        error('test_salpha_Indented_Miller_kappascan:BadFluxDerivative', ...
            'The signed dpsi/da must be finite and nonzero.');
    end
    Vprime = dVda/merluc.dpsi_da;

    if any(~isfinite([V, Vprime, dVda])) ...
            || V <= 0 || dVda <= 0 || Vprime == 0
        error('test_salpha_Indented_Miller_kappascan:BadVolumeIntegral', ...
            'A modified V, Vprime, or dV/da value is invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
