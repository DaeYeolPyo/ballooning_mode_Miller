clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Reference equilibrium and indentation scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_I0.8');
target_psiN = 0.7;
I0_grid = 0:0.2:0.6;

ntheta0 = 16;
nshat = 30;
shat_end = 7.0;
nalpha = 30;
alpha_end = 12.0;
theta_bnd = 5*pi;
n_surface_points = int32(1024);
n_fit_points = int32(256);
n_balloon_points = 151;

%% Fit the reference surface once
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

nI0 = numel(I0_grid);
lambda_max_I0 = nan(nI0, nshat, nalpha);
best_theta0_I0 = nan(nI0, nshat, nalpha);
lambda_theta0_I0 = nan(nI0, nshat, nalpha, ntheta0);
fail_message = strings(nI0, nshat, nalpha, ntheta0);

volume_I0 = nan(1, nI0);
Vprime_I0 = nan(1, nI0);
dVda_I0 = nan(1, nI0);
dpsi_da_I0 = nan(1, nI0);
q_check_I0 = nan(1, nI0);
minimum_abs_jacobian = nan(1, nI0);
qprime_scan = nan(nI0, nshat);
pprime_scan = nan(nI0, nalpha);
shat_equilibrium = nan(1, nI0);
alpha_equilibrium = nan(1, nI0);
param_I0 = cell(1, nI0);
merluc_I0 = cell(1, nI0);

fprintf('Indented-Miller indentation scan reference\n');
fprintf('  GEQDSK                      : %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  reference I0                : %.8f\n', ...
    param_reference.I0);
fprintf('  reference s_I0=a*dI0/da     : %+.8f\n', ...
    param_reference.s_I0);
fprintf('  reference delta / kappa     : %+.8f / %.8f\n', ...
    param_reference.delta, param_reference.kappa);
fprintf('  reference Bp relative RMS   : %.6e\n', ...
    Bp_reference_relative_rms);
fprintf('  scan I0 values              :');
fprintf(' %.1f', I0_grid);
fprintf('\n');

assert(Bp_reference_relative_rms < 2.e-2, ...
    'The reference indented-Miller Bp fit exceeds 2%% RMS error.');

%% Build all modified local equilibria
for ii = 1:nI0
    I0 = I0_grid(ii);
    param_i = set_indentation( ...
        param_reference, I0, param_reference.s_I0);
    bnd_i = Indented_Miller_boundary( ...
        param_i.R0, param_i.a, param_i.delta, ...
        param_i.kappa, param_i.I0, NTheta=n_fit_points);
    merluc_i = Indented_Miller_Mercier_Luc(surf, param_i, bnd_i);

    [volume_I0(ii), Vprime_I0(ii), dVda_I0(ii)] = ...
        indented_volume_integrals(merluc_i);
    dpsi_da_I0(ii) = merluc_i.dpsi_da;
    q_check_I0(ii) = merluc_i.q_check;
    minimum_abs_jacobian(ii) = min(abs(merluc_i.jac));

    normal_projection_a = ...
        merluc_i.normal_sign*merluc_i.normal_projection_raw;
    normal_projection_psi = ...
        sign(merluc_i.dpsi_da)*normal_projection_a;
    if any(normal_projection_psi <= 0)
        error('test_salpha_Indented_Miller_Iscan:BadNormal', ...
            'The Mercier normal is misaligned at I0=%g.', I0);
    end

    jacobian_scale = max(abs(merluc_i.jac));
    if minimum_abs_jacobian(ii) <= 1.e-10*jacobian_scale
        error('test_salpha_Indented_Miller_Iscan:SingularJacobian', ...
            'The fitted coordinate is nearly singular at I0=%g.', I0);
    end
    if abs(q_check_I0(ii) - surf.q) ...
            > 1.e-10*max(1, abs(surf.q))
        error('test_salpha_Indented_Miller_Iscan:BadQ', ...
            'The modified surface does not reproduce q at I0=%g.', I0);
    end

    qprime_scan(ii, :) = shat.*( ...
        surf.q*Vprime_I0(ii)/(2*volume_I0(ii)));
    pprime_scan(ii, :) = -alpha.*(4*pi^2/(2*Vprime_I0(ii))) ...
        .*sqrt(2*pi^2*eq.rmaxis/volume_I0(ii))/mu0;
    shat_equilibrium(ii) = 2*volume_I0(ii)*qprime_equilibrium ...
        /(surf.q*Vprime_I0(ii));
    alpha_equilibrium(ii) = -mu0*pprime_equilibrium ...
        *(2*Vprime_I0(ii))/(4*pi^2) ...
        *sqrt(volume_I0(ii)/(2*pi^2*eq.rmaxis));

    param_I0{ii} = param_i;
    merluc_I0{ii} = merluc_i;

    fprintf(['  I0=%.1f: V=%.6e, Vprime=%+.6e, ', ...
        'dpsi/da=%+.6e, min|J|=%.3e\n'], ...
        I0, volume_I0(ii), Vprime_I0(ii), ...
        dpsi_da_I0(ii), minimum_abs_jacobian(ii));
end

%% I0 -> shat -> alpha -> theta0 scan
nsolve = nI0*nshat*nalpha*ntheta0;
fprintf('Reference setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting indentation scan: %d I0 x %d shat x %d alpha x ', ...
    '%d theta0 = %d solves.\n'], ...
    nI0, nshat, nalpha, ntheta0, nsolve);

scan_tic = tic;
total_rows = nI0*nshat;

for ii = 1:nI0
    indentation_tic = tic;
    merluc_i = merluc_I0{ii};
    fprintf('\nI0 %d/%d: %.1f\n', ii, nI0, I0_grid(ii));

    for ss = 1:nshat
        row_tic = tic;
        qprime = qprime_scan(ii, ss);

        for aa = 1:nalpha
            pprime = pprime_scan(ii, aa);
            [~, ~, theta_diag] = evaluate_theta_psi(merluc_i, ...
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
                    lambda_theta0_I0(ii, ss, aa, tt) = lambda(1);
                catch ME
                    fail_message(ii, ss, aa, tt) = string(ME.message);
                end
            end

            values = reshape( ...
                lambda_theta0_I0(ii, ss, aa, :), 1, []);
            [lambda_max_I0(ii, ss, aa), imax] = ...
                max(values, [], 'omitnan');
            if isfinite(lambda_max_I0(ii, ss, aa))
                best_theta0_I0(ii, ss, aa) = theta0(imax);
            end
        end

        completed_rows = (ii - 1)*nshat + ss;
        elapsed = toc(scan_tic);
        eta = elapsed/completed_rows*(total_rows - completed_rows);
        fprintf(['  shat row %2d/%2d: shat=%.4g, row %.1f s, ', ...
            'total ETA %.1f min\n'], ss, nshat, shat(ss), ...
            toc(row_tic), eta/60);
    end

    fprintf('I0=%.1f finished in %.1f min.\n', ...
        I0_grid(ii), toc(indentation_tic)/60);
end

elapsed = toc(scan_tic);
nfailed = nnz(strlength(fail_message) > 0);
fprintf(['\nAll indentation scans finished in %.1f min ', ...
    '(%d failed solves).\n'], elapsed/60, nfailed);

%% Package and validate results
salpha_Indented_Miller_Iscan_result = struct();
salpha_Indented_Miller_Iscan_result.model = ...
    'indented-miller-I0-scan-mercier-luc';
salpha_Indented_Miller_Iscan_result.I0 = I0_grid;
salpha_Indented_Miller_Iscan_result.s_I0 = ...
    param_reference.s_I0;
salpha_Indented_Miller_Iscan_result.shat = shat;
salpha_Indented_Miller_Iscan_result.alpha = alpha;
salpha_Indented_Miller_Iscan_result.theta0 = theta0;
salpha_Indented_Miller_Iscan_result.lambda_max = lambda_max_I0;
salpha_Indented_Miller_Iscan_result.lambda_theta0 = ...
    lambda_theta0_I0;
salpha_Indented_Miller_Iscan_result.best_theta0 = best_theta0_I0;
salpha_Indented_Miller_Iscan_result.fail_message = fail_message;
salpha_Indented_Miller_Iscan_result.qprime = qprime_scan;
salpha_Indented_Miller_Iscan_result.pprime = pprime_scan;
salpha_Indented_Miller_Iscan_result.volume = volume_I0;
salpha_Indented_Miller_Iscan_result.Vprime = Vprime_I0;
salpha_Indented_Miller_Iscan_result.dVda = dVda_I0;
salpha_Indented_Miller_Iscan_result.dpsi_da = dpsi_da_I0;
salpha_Indented_Miller_Iscan_result.q_check = q_check_I0;
salpha_Indented_Miller_Iscan_result.minimum_abs_jacobian = ...
    minimum_abs_jacobian;
salpha_Indented_Miller_Iscan_result.shat_equilibrium = ...
    shat_equilibrium;
salpha_Indented_Miller_Iscan_result.alpha_equilibrium = ...
    alpha_equilibrium;
salpha_Indented_Miller_Iscan_result.reference_param = ...
    param_reference;
salpha_Indented_Miller_Iscan_result.param = param_I0;
salpha_Indented_Miller_Iscan_result.target_psiN = target_psiN;
salpha_Indented_Miller_Iscan_result.geqdsk_file = geqdsk_file;
salpha_Indented_Miller_Iscan_result.Bp_reference_relative_rms = ...
    Bp_reference_relative_rms;
salpha_Indented_Miller_Iscan_result.scan_elapsed_seconds = elapsed;
salpha_Indented_Miller_Iscan_result.nfailed = nfailed;

has_marginal_crossing = false(1, nI0);
lambda_range = nan(nI0, 2);
for ii = 1:nI0
    lambda_i = squeeze(lambda_max_I0(ii, :, :));
    finite_i = lambda_i(isfinite(lambda_i));
    if ~isempty(finite_i)
        lambda_range(ii, :) = [min(finite_i), max(finite_i)];
        has_marginal_crossing(ii) = lambda_range(ii, 1) <= 0 ...
            && lambda_range(ii, 2) >= 0;
    end
end
salpha_Indented_Miller_Iscan_result.lambda_range = lambda_range;
salpha_Indented_Miller_Iscan_result.has_marginal_crossing = ...
    has_marginal_crossing;

assert(nfailed == 0, ...
    'The indentation scan contains failed eigenvalue solves.');
assert(all(isfinite(lambda_max_I0), 'all'), ...
    'The indentation scan contains non-finite eigenvalues.');
assert(all(lambda_max_I0(:, :, 1) < 0, 'all'), ...
    'The alpha=0 column should be stable for every I0.');

fprintf('\nIndentation scan summary\n');
for ii = 1:nI0
    fprintf(['  I0=%.1f: lambda=[%+.4e,%+.4e], crossing=%d, ', ...
        '(shat,alpha)_eq=(%.4f,%.4f)\n'], ...
        I0_grid(ii), lambda_range(ii, 1), lambda_range(ii, 2), ...
        has_marginal_crossing(ii), ...
        shat_equilibrium(ii), alpha_equilibrium(ii));
end

%% Individual stability maps with a shared color scale
finite_lambda = lambda_max_I0(isfinite(lambda_max_I0));
color_limits = [min(finite_lambda), max(finite_lambda)];
map_figure = figure('Color', 'w', 'Name', ...
    'Indented-Miller I0 scan stability maps');
map_layout = tiledlayout(map_figure, 2, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
map_axes = gobjects(1, nI0);

for ii = 1:nI0
    map_axes(ii) = nexttile(map_layout);
    lambda_i = squeeze(lambda_max_I0(ii, :, :));
    contourf(map_axes(ii), alpha, shat, lambda_i, 30, ...
        'LineStyle', 'none');
    hold(map_axes(ii), 'on');
    if has_marginal_crossing(ii)
        contour(map_axes(ii), alpha, shat, lambda_i, [0, 0], ...
            'k', 'LineWidth', 1.8);
    end
    plot(map_axes(ii), alpha_equilibrium(ii), ...
        shat_equilibrium(ii), 'p', ...
        'MarkerSize', 9, 'MarkerEdgeColor', 'w', ...
        'MarkerFaceColor', [0.85, 0.1, 0.75], 'LineWidth', 1.0);
    clim(map_axes(ii), color_limits);
    grid(map_axes(ii), 'on');
    box(map_axes(ii), 'on');
    map_axes(ii).Layer = 'top';
    xlabel(map_axes(ii), '$\alpha$', 'Interpreter', 'latex');
    ylabel(map_axes(ii), '$\hat{s}$', 'Interpreter', 'latex');
    title(map_axes(ii), sprintf('$I_0=%.1f$', I0_grid(ii)), ...
        'Interpreter', 'latex');
end

cb_map = colorbar(map_axes(end));
cb_map.Label.Interpreter = 'latex';
cb_map.Label.String = '$\lambda_{\max}$';

%% Overlay marginal-stability contours
marginal_figure = figure('Color', 'w', 'Name', ...
    'Indented-Miller I0 marginal curves');
ax_marginal = axes(marginal_figure);
hold(ax_marginal, 'on');
colors = lines(nI0);
curve_handles = gobjects(0);
curve_labels = strings(0);

for ii = 1:nI0
    lambda_i = squeeze(lambda_max_I0(ii, :, :));
    if has_marginal_crossing(ii)
        [~, curve] = contour(ax_marginal, alpha, shat, ...
            lambda_i, [0, 0], 'Color', colors(ii, :), ...
            'LineWidth', 2.2);
        curve_handles(end + 1) = curve; %#ok<SAGROW>
        curve_labels(end + 1) = ...
            sprintf('$I_0=%.1f$', I0_grid(ii)); %#ok<SAGROW>
    else
        warning('test_salpha_Indented_Miller_Iscan:NoCrossing', ...
            ['No marginal crossing was found for I0=%g. ', ...
             'Increase the scan range.'], I0_grid(ii));
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

salpha_Indented_Miller_Iscan_result.map_figure = map_figure;
salpha_Indented_Miller_Iscan_result.marginal_figure = ...
    marginal_figure;
salpha_Indented_Miller_Iscan_result.status = 'PASS';

fprintf('  validation status           : PASS\n');

function param = set_indentation(param, I0, s_I0)
    param.I0 = I0;
    param.indentation = I0;
    param.d = param.a*I0;
    param.coeff(5) = param.d;

    % Preserve the reference dimensionless radial indentation shear
    % s_I0=a*dI0/da while changing the local value of I0.
    param.s_I0 = s_I0;
    param.s_I = s_I0;
    param.dI0_da = s_I0/param.a;
    param.dI0_dr = param.dI0_da;
    param.dd_da = I0 + s_I0;
    param.dcoeff_da(5) = param.dd_da;
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
        error('test_salpha_Indented_Miller_Iscan:BadFluxDerivative', ...
            'The signed dpsi/da must be finite and nonzero.');
    end
    Vprime = dVda/merluc.dpsi_da;

    if any(~isfinite([V, Vprime, dVda])) ...
            || V <= 0 || dVda <= 0 || Vprime == 0
        error('test_salpha_Indented_Miller_Iscan:BadVolumeIntegral', ...
            'A modified V, Vprime, or dV/da value is invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
