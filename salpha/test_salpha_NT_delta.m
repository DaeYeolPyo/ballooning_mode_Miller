clc
close all
clear

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% NT reference equilibrium and delta-scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_NT0.6');

if ~exist('target_psiN', 'var') || isempty(target_psiN)
    target_psiN = 0.8;
end
if ~exist('miller_fit_dpsi', 'var') || isempty(miller_fit_dpsi)
    miller_fit_dpsi = 8.e-2;
end
if ~exist('delta_grid', 'var') || isempty(delta_grid)
    delta_grid = [-0.6, -0.4, -0.2];
end
if ~exist('sdelta_mode', 'var') || isempty(sdelta_mode)
    sdelta_mode = "reference";
end
if ~exist('ntheta0', 'var') || isempty(ntheta0)
    ntheta0 = 16;
end
if ~exist('nshat', 'var') || isempty(nshat)
    nshat = 30;
end
if ~exist('shat_end', 'var') || isempty(shat_end)
    shat_end = 7.0;
end
if ~exist('nalpha', 'var') || isempty(nalpha)
    nalpha = 30;
end
if ~exist('alpha_end', 'var') || isempty(alpha_end)
    alpha_end = 12.0;
end
if ~exist('theta_bnd', 'var') || isempty(theta_bnd)
    theta_bnd = 5*pi;
end
if ~exist('n_surface_points', 'var') || isempty(n_surface_points)
    n_surface_points = 1024;
end
if ~exist('n_balloon_points', 'var') || isempty(n_balloon_points)
    n_balloon_points = 151;
end

delta_grid = delta_grid(:).';
if any(~isfinite(delta_grid)) || any(abs(delta_grid) >= 1)
    error('test_salpha_NT_delta:BadDelta', ...
        'Every delta value must be finite and satisfy abs(delta) < 1.');
end
sdelta_mode = lower(string(sdelta_mode));
if ~ismember(sdelta_mode, ["reference", "zero"])
    error('test_salpha_NT_delta:BadSdeltaMode', ...
        'sdelta_mode must be "reference" or "zero".');
end

%% Fit the reference NT surface once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=int32(n_surface_points));
[param_reference, bnd_reference] = fit_Miller( ...
    eq, eqfunc, target_psiN, NTheta=int32(n_surface_points), ...
    dpsi=miller_fit_dpsi);
merluc_reference = Miller_Mercier_Luc( ...
    surf, param_reference, bnd_reference);

eq_on_reference = eqfunc.eval( ...
    merluc_reference.R, merluc_reference.Z);
Bp_reference_relative_rms = relative_rms( ...
    merluc_reference.Bp, eq_on_reference.Bp);

FFprime = eqfunc.FFprime(target_psiN);
theta0 = 2*pi*((0:ntheta0-1) - floor(ntheta0/2))/ntheta0;
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
sdelta_used = nan(1, ndelta);
qprime_scan = nan(ndelta, nshat);
pprime_scan = nan(ndelta, nalpha);
param_delta = cell(1, ndelta);

fprintf('NT delta-scan Miller reference\n');
fprintf('  GEQDSK                    : %s\n', geqdsk_file);
fprintf('  target psiN               : %.6f\n', target_psiN);
fprintf('  reference delta/kappa     : %+.8f / %.8f\n', ...
    param_reference.delta, param_reference.kappa);
fprintf('  reference s_delta         : %+.8f\n', ...
    param_reference.s_delta);
fprintf('  radial fit half width     : %.6f in psiN\n', ...
    param_reference.radial_fit.half_width);
fprintf('  reference Bp relative RMS : %.6e\n', ...
    Bp_reference_relative_rms);
fprintf('  s_delta scan mode         : %s\n', sdelta_mode);
fprintf('  delta values              :');
fprintf(' %+.4f', delta_grid);
fprintf('\n');

assert(param_reference.delta < 0, ...
    'The fitted reference equilibrium is not negative triangularity.');
assert(Bp_reference_relative_rms < 5.e-2, ...
    'The reference NT Miller Bp fit exceeds 5%% RMS error.');

%% delta -> shat -> alpha -> theta0 scan
nsolve = ndelta*nshat*nalpha*ntheta0;
fprintf('Reference setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting NT delta scan: %d delta x %d shat x %d alpha x ', ...
    '%d theta0 = %d solves.\n'], ...
    ndelta, nshat, nalpha, ntheta0, nsolve);

scan_tic = tic;
total_rows = ndelta*nshat;

for dd = 1:ndelta
    delta_tic = tic;
    param_d = param_reference;
    param_d.delta = delta_grid(dd);
    if sdelta_mode == "zero"
        param_d.s_delta = 0.0;
    end
    sdelta_used(dd) = param_d.s_delta;

    bnd_d = Miller_boundary(param_d.R0, param_d.r, ...
        param_d.delta, param_d.kappa, ...
        NTheta=int32(n_surface_points));
    merluc_d = Miller_Mercier_Luc(surf, param_d, bnd_d);
    param_delta{dd} = param_d;

    [volume_delta(dd), Vprime_delta(dd), dVdr_delta(dd)] = ...
        miller_volume_integrals(merluc_d);
    q_check_delta(dd) = merluc_d.q_check;
    minimum_abs_jacobian(dd) = min(abs(merluc_d.jac));

    jacobian_scale = max(abs(merluc_d.jac));
    if minimum_abs_jacobian(dd) <= 1.e-10*jacobian_scale
        error('test_salpha_NT_delta:NearlySingularJacobian', ...
            'The Miller Jacobian is nearly singular at delta=%g.', ...
            delta_grid(dd));
    end
    if abs(q_check_delta(dd) - surf.q) ...
            > 1.e-10*max(1, abs(surf.q))
        error('test_salpha_NT_delta:BadQ', ...
            'The modified Miller surface does not reproduce q.');
    end

    qprime_scan(dd, :) = shat.*( ...
        surf.q*Vprime_delta(dd)/(2*volume_delta(dd)));
    pprime_scan(dd, :) = -alpha.*(4*pi^2/(2*Vprime_delta(dd))) ...
        .*sqrt(2*pi^2*eq.rmaxis/volume_delta(dd))/mu0;

    fprintf(['\nDelta %d/%d: delta=%+.4f, s_delta=%+.4f, ', ...
        'V=%.6g, Vprime=%+.6g, min|J|=%.3e\n'], ...
        dd, ndelta, delta_grid(dd), sdelta_used(dd), ...
        volume_delta(dd), Vprime_delta(dd), ...
        minimum_abs_jacobian(dd));

    for ss = 1:nshat
        row_tic = tic;
        qprime = qprime_scan(dd, ss);

        for aa = 1:nalpha
            pprime = pprime_scan(dd, aa);
            [~, ~, theta_diag] = evaluate_theta_psi(merluc_d, ...
                FFprime, pprime, qprime, EnforcePeriodicity=true);

            for tt = 1:ntheta0
                try
                    bal = evaluate_gcf(theta_diag, theta0(tt), ...
                        theta_bnd=theta_bnd, a_N=aN, B_N=BN);
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
        fprintf(['  shat row %2d/%2d: shat=%.4g, row %.1f s, ', ...
            'total ETA %.1f min\n'], ss, nshat, shat(ss), ...
            toc(row_tic), eta/60);
    end

    fprintf('Delta %+.4f finished in %.1f min.\n', ...
        delta_grid(dd), toc(delta_tic)/60);
end

elapsed = toc(scan_tic);
nfailed = nnz(strlength(fail_message) > 0);
fprintf('\nNT delta scans finished in %.1f min (%d failed solves).\n', ...
    elapsed/60, nfailed);

%% Package results and validate
salpha_NT_delta_result = struct();
salpha_NT_delta_result.model = ...
    'negative-triangularity-delta-scan-miller-mercier-luc';
salpha_NT_delta_result.delta = delta_grid;
salpha_NT_delta_result.sdelta_mode = sdelta_mode;
salpha_NT_delta_result.sdelta_used = sdelta_used;
salpha_NT_delta_result.shat = shat;
salpha_NT_delta_result.alpha = alpha;
salpha_NT_delta_result.theta0 = theta0;
salpha_NT_delta_result.lambda_max = lambda_max_delta;
salpha_NT_delta_result.lambda_theta0 = lambda_theta0_delta;
salpha_NT_delta_result.best_theta0 = best_theta0_delta;
salpha_NT_delta_result.fail_message = fail_message;
salpha_NT_delta_result.qprime = qprime_scan;
salpha_NT_delta_result.pprime = pprime_scan;
salpha_NT_delta_result.volume = volume_delta;
salpha_NT_delta_result.Vprime = Vprime_delta;
salpha_NT_delta_result.dVdr = dVdr_delta;
salpha_NT_delta_result.q_check = q_check_delta;
salpha_NT_delta_result.minimum_abs_jacobian = minimum_abs_jacobian;
salpha_NT_delta_result.reference_param = param_reference;
salpha_NT_delta_result.param = param_delta;
salpha_NT_delta_result.target_psiN = target_psiN;
salpha_NT_delta_result.geqdsk_file = geqdsk_file;
salpha_NT_delta_result.Bp_reference_relative_rms = ...
    Bp_reference_relative_rms;
salpha_NT_delta_result.scan_elapsed_seconds = elapsed;
salpha_NT_delta_result.nfailed = nfailed;

has_crossing = false(1, ndelta);
for dd = 1:ndelta
    lambda_d = squeeze(lambda_max_delta(dd, :, :));
    has_crossing(dd) = min(lambda_d, [], 'all', 'omitnan') <= 0 ...
        && max(lambda_d, [], 'all', 'omitnan') >= 0;
end
salpha_NT_delta_result.has_marginal_crossing = has_crossing;

assert(nfailed == 0, ...
    'The NT delta scan contains failed eigenvalue solves.');
assert(all(isfinite(lambda_max_delta), 'all'), ...
    'The NT delta scan contains non-finite eigenvalues.');
assert(all(lambda_max_delta(:, :, 1) < 0, 'all'), ...
    'The alpha=0 column should be stable for every delta.');

fprintf('  marginal crossing by delta :');
fprintf(' %d', has_crossing);
fprintf('\n');

%% Overlay marginal-stability curves
figure('Color', 'w', 'Name', 'NT delta s-alpha marginal curves');
ax = axes();
hold(ax, 'on');
colors = lines(ndelta);
curve_handles = gobjects(0);
curve_labels = strings(0);

for dd = 1:ndelta
    lambda_d = squeeze(lambda_max_delta(dd, :, :));
    if has_crossing(dd)
        [~, h] = contour(ax, alpha, shat, lambda_d, [0, 0], ...
            'Color', colors(dd, :), 'LineWidth', 2.2);
        curve_handles(end + 1) = h; %#ok<SAGROW>
        curve_labels(end + 1) = ...
            sprintf('$\\delta=%+.3g$', delta_grid(dd)); %#ok<SAGROW>
    else
        warning('test_salpha_NT_delta:NoMarginalCrossing', ...
            ['No lambda_max=0 crossing for delta=%g. Increase ', ...
             'shat_end or alpha_end.'], delta_grid(dd));
    end
end

grid(ax, 'on');
box(ax, 'on');
ax.Layer = 'top';
ax.TickDir = 'out';
xlabel(ax, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax, '$\hat{s}$', 'Interpreter', 'latex');
if ~isempty(curve_handles)
    legend(ax, curve_handles, curve_labels, ...
        'Interpreter', 'latex', 'Location', 'best');
end

fprintf('  validation status          : PASS\n');

function [V, Vprime, dVdr] = miller_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVdr = 2*pi*trapz(theta, R.*abs(jac));
    if ~isfinite(merluc.dpsi_dr) || abs(merluc.dpsi_dr) <= eps
        error('test_salpha_NT_delta:BadFluxDerivative', ...
            'The Miller dpsi/dr must be finite and nonzero.');
    end
    Vprime = dVdr/merluc.dpsi_dr;
    if any(~isfinite([V, Vprime, dVdr])) ...
            || V <= 0 || dVdr <= 0 || Vprime == 0
        error('test_salpha_NT_delta:BadVolumeIntegral', ...
            'The Miller volume derivatives are invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
