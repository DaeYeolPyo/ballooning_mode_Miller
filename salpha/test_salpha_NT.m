clc
close all
clear

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Negative-triangularity equilibrium and scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_NT0.6');

if ~exist('target_psiN', 'var') || isempty(target_psiN)
    target_psiN = 0.7;
end
if ~exist('miller_fit_dpsi', 'var') || isempty(miller_fit_dpsi)
    miller_fit_dpsi = 8.e-2;
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

validateattributes(target_psiN, {'double'}, ...
    {'scalar', 'finite', '>', 0, '<', 1});
validateattributes(miller_fit_dpsi, {'double'}, ...
    {'scalar', 'finite', 'positive'});
validateattributes(ntheta0, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1});
validateattributes(nshat, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2});
validateattributes(nalpha, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2});
validateattributes(n_surface_points, {'numeric'}, ...
    {'scalar', 'integer', '>=', 64});
validateattributes(n_balloon_points, {'numeric'}, ...
    {'scalar', 'integer', '>=', 21});

%% Fit the NT Miller local equilibrium once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=int32(n_surface_points));
ints = compute_contour_integrals(eq, eqfunc, surf, ...
    isNormalized=false);

[param, bnd] = fit_Miller(eq, eqfunc, target_psiN, ...
    NTheta=int32(n_surface_points), dpsi=miller_fit_dpsi);
merluc = Miller_Mercier_Luc(surf, param, bnd);
[volume, Vprime, dVdr] = miller_volume_integrals(merluc);

normal_projection_r = ...
    merluc.normal_sign*merluc.normal_projection_raw;
normal_projection_psi = sign(merluc.dpsi_dr)*normal_projection_r;
assert(all(normal_projection_psi > 0), ...
    'The NT Miller normal is not aligned with grad(psi).');

eq_on_miller = eqfunc.eval(merluc.R, merluc.Z);
Bp_fit_relative_rms = relative_rms(merluc.Bp, eq_on_miller.Bp);
psiN_fit_rms = sqrt(mean((eq_on_miller.psiN(:) - target_psiN).^2));
volume_relative_error = abs(volume - ints.V)/max(abs(ints.V), eps);
Vprime_relative_error = abs(abs(Vprime) - abs(ints.Vprime)) ...
    /max(abs(ints.Vprime), eps);

FFprime = eqfunc.FFprime(target_psiN);
theta0 = 2*pi*((0:ntheta0-1) - floor(ntheta0/2))/ntheta0;
alpha = linspace(0, alpha_end, nalpha);
shat = linspace(0, shat_end, nshat);
theta_grid = linspace(-theta_bnd, theta_bnd, n_balloon_points).';

aN = eq.rmaxis;
BN = abs(eq.bcentr);
mu0 = 4*pi*1.e-7;

fprintf('Negative-triangularity Miller s-alpha setup\n');
fprintf('  GEQDSK                      : %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  R0, r, kappa, delta        : %.8e, %.8e, %.8e, %+.8e\n', ...
    param.R0, param.r, param.kappa, param.delta);
fprintf('  dR0/dr, s_kappa, s_delta   : %+.8e, %+.8e, %+.8e\n', ...
    param.dR0_dr, param.s_kappa, param.s_delta);
fprintf('  radial fit half width      : %.6f in psiN\n', ...
    param.radial_fit.half_width);
fprintf('  Miller-path psiN RMS       : %.6e\n', psiN_fit_rms);
fprintf('  Bp relative RMS            : %.6e\n', Bp_fit_relative_rms);
fprintf('  q target / check           : %.10f / %.10f\n', ...
    surf.q, merluc.q_check);
fprintf('  V Miller / GEQDSK          : %.8e / %.8e m^3\n', ...
    volume, ints.V);
fprintf('  |Vprime| Miller / GEQDSK   : %.8e / %.8e m^3/Wb\n', ...
    abs(Vprime), abs(ints.Vprime));
fprintf('  V / Vprime relative errors : %.6e / %.6e\n', ...
    volume_relative_error, Vprime_relative_error);
fprintf('  dV/dr, signed dpsi/dr      : %.8e, %+.8e\n', ...
    dVdr, merluc.dpsi_dr);

assert(param.delta < 0, ...
    'The fitted surface is not negative triangularity.');
assert(Bp_fit_relative_rms < 5.e-2, ...
    'The NT Miller Bp fit exceeds the 5%% RMS validation threshold.');
assert(abs(merluc.q_check - surf.q) < 1.e-10*max(1, abs(surf.q)), ...
    'The NT Miller local equilibrium does not reproduce q.');
assert(volume_relative_error < 5.e-2 ...
        && Vprime_relative_error < 5.e-2, ...
    'The NT Miller volume normalization differs from GEQDSK by over 5%%.');

%% shat -> alpha -> theta0 scan
lambda_max = nan(nshat, nalpha);
best_theta0 = nan(nshat, nalpha);
lambda_theta0 = nan(nshat, nalpha, ntheta0);
fail_message = strings(nshat, nalpha, ntheta0);
qprime_scan = shat.*(surf.q*Vprime/(2*volume));
pprime_scan = -alpha.*(4*pi^2/(2*Vprime)) ...
    .*sqrt(2*pi^2*eq.rmaxis/volume)/mu0;

nsolve = nshat*nalpha*ntheta0;
fprintf('Equilibrium setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting NT s-alpha scan: %d shat x %d alpha x ', ...
    '%d theta0 = %d solves.\n'], nshat, nalpha, ntheta0, nsolve);

scan_tic = tic;
for ss = 1:nshat
    row_tic = tic;
    qprime = qprime_scan(ss);

    for aa = 1:nalpha
        pprime = pprime_scan(aa);
        [~, ~, theta_diag] = evaluate_theta_psi(merluc, ...
            FFprime, pprime, qprime, EnforcePeriodicity=true);

        for tt = 1:ntheta0
            try
                bal = evaluate_gcf(theta_diag, theta0(tt), ...
                    theta_bnd=theta_bnd, a_N=aN, B_N=BN);

                % Preserve the non-periodic shear term on the reduced
                % finite-element grid used for the repeated solves.
                g = interp1(bal.theta, bal.g, theta_grid, 'pchip');
                c = interp1(bal.theta, bal.c, theta_grid, 'pchip');
                f = interp1(bal.theta, bal.f, theta_grid, 'pchip');

                [Gmat, Cmat, Fmat, theta_Dof] = construct_matrix( ...
                    theta_grid, g, c, f);
                [lambda, ~] = calculate_eigenvalue( ...
                    theta_Dof, Gmat, Cmat, Fmat, 'dirichlet');
                lambda_theta0(ss, aa, tt) = lambda(1);
            catch ME
                fail_message(ss, aa, tt) = string(ME.message);
            end
        end

        values = reshape(lambda_theta0(ss, aa, :), 1, []);
        [lambda_max(ss, aa), imax] = max(values, [], 'omitnan');
        if isfinite(lambda_max(ss, aa))
            best_theta0(ss, aa) = theta0(imax);
        end
    end

    elapsed = toc(scan_tic);
    eta = elapsed/ss*(nshat - ss);
    fprintf(['  shat row %2d/%2d: shat = %.4g, row %.1f s, ', ...
        'ETA %.1f min\n'], ss, nshat, shat(ss), ...
        toc(row_tic), eta/60);
end

elapsed = toc(scan_tic);
nfailed = nnz(strlength(fail_message) > 0);
fprintf('NT s-alpha scan finished in %.1f min (%d failed solves).\n', ...
    elapsed/60, nfailed);

%% Package scan products
salpha_NT_result = struct();
salpha_NT_result.model = 'negative-triangularity-miller-mercier-luc';
salpha_NT_result.shat = shat;
salpha_NT_result.alpha = alpha;
salpha_NT_result.theta0 = theta0;
salpha_NT_result.lambda_max = lambda_max;
salpha_NT_result.lambda_theta0 = lambda_theta0;
salpha_NT_result.best_theta0 = best_theta0;
salpha_NT_result.fail_message = fail_message;
salpha_NT_result.qprime = qprime_scan;
salpha_NT_result.pprime = pprime_scan;
salpha_NT_result.target_psiN = target_psiN;
salpha_NT_result.geqdsk_file = geqdsk_file;
salpha_NT_result.param = param;
salpha_NT_result.volume = volume;
salpha_NT_result.Vprime = Vprime;
salpha_NT_result.dVdr = dVdr;
salpha_NT_result.volume_relative_error = volume_relative_error;
salpha_NT_result.Vprime_relative_error = Vprime_relative_error;
salpha_NT_result.Bp_fit_relative_rms = Bp_fit_relative_rms;
salpha_NT_result.psiN_fit_rms = psiN_fit_rms;
salpha_NT_result.scan_elapsed_seconds = elapsed;
salpha_NT_result.nfailed = nfailed;

finite_lambda = lambda_max(isfinite(lambda_max));
if isempty(finite_lambda)
    error('test_salpha_NT:NoFiniteEigenvalues', ...
        'The NT s-alpha scan produced no finite eigenvalues.');
end
has_marginal_crossing = min(finite_lambda) <= 0 ...
    && max(finite_lambda) >= 0;
salpha_NT_result.has_marginal_crossing = has_marginal_crossing;

fprintf('  lambda_max range            : [%+.8e, %+.8e]\n', ...
    min(finite_lambda), max(finite_lambda));
fprintf('  marginal crossing found     : %d\n', has_marginal_crossing);

assert(nfailed == 0, ...
    'The NT s-alpha scan contains failed eigenvalue solves.');
assert(all(isfinite(lambda_max), 'all'), ...
    'The NT s-alpha map contains non-finite eigenvalues.');
assert(all(lambda_max(:, 1) < 0), ...
    'The alpha=0 column should be stable on the finite Dirichlet domain.');

%% Stability and marginal-contour plot
figure('Color', 'w', 'Name', 'NT s-alpha stability diagram');
ax_stability = axes();
contourf(ax_stability, alpha, shat, lambda_max, 30, ...
    'LineStyle', 'none');
hold(ax_stability, 'on');
if has_marginal_crossing
    contour(ax_stability, alpha, shat, lambda_max, [0, 0], ...
        'k', 'LineWidth', 2.2);
else
    warning('test_salpha_NT:NoMarginalCrossing', ...
        ['No lambda_max=0 crossing was found. Increase shat_end or ', ...
         'alpha_end, or refine the scan range.']);
end
cb_stability = colorbar(ax_stability);
cb_stability.Label.Interpreter = 'latex';
cb_stability.Label.String = '$\lambda_{\max}$';
grid(ax_stability, 'on');
box(ax_stability, 'on');
ax_stability.Layer = 'top';
xlabel(ax_stability, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax_stability, '$\hat{s}$', 'Interpreter', 'latex');

%% Ballooning phase that maximizes lambda
[alpha_grid, shat_grid] = meshgrid(alpha, shat);
valid = isfinite(best_theta0);

figure('Color', 'w', 'Name', ...
    'NT phase of maximum ballooning eigenvalue');
ax_phase = axes();
scatter(ax_phase, alpha_grid(valid), shat_grid(valid), 58, ...
    best_theta0(valid)/pi, 's', 'filled', ...
    'MarkerEdgeColor', [0.18, 0.18, 0.18], 'LineWidth', 0.25);

theta0_over_pi = theta0/pi;
dtheta0_over_pi = 2/ntheta0;
colormap(ax_phase, hsv(ntheta0));
clim(ax_phase, [theta0_over_pi(1) - dtheta0_over_pi/2, ...
                theta0_over_pi(end) + dtheta0_over_pi/2]);
cb_phase = colorbar(ax_phase);
cb_phase.Ticks = theta0_over_pi;
cb_phase.TickLabels = compose('%+.2f', theta0_over_pi);
cb_phase.Label.Interpreter = 'latex';
cb_phase.Label.String = '$\theta_{0,\mathrm{max}}/\pi$';

dalpha = alpha(2) - alpha(1);
dshat = shat(2) - shat(1);
xlim(ax_phase, [min(alpha) - 0.5*dalpha, ...
                max(alpha) + 0.5*dalpha]);
ylim(ax_phase, [min(shat) - 0.5*dshat, ...
                max(shat) + 0.5*dshat]);
grid(ax_phase, 'on');
box(ax_phase, 'on');
ax_phase.Layer = 'top';
ax_phase.TickDir = 'out';
xlabel(ax_phase, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax_phase, '$\hat{s}$', 'Interpreter', 'latex');
title(ax_phase, 'NT phase of maximum ballooning eigenvalue');

fprintf('  validation status           : PASS\n');

function [V, Vprime, dVdr] = miller_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVdr = 2*pi*trapz(theta, R.*abs(jac));

    if ~isfinite(merluc.dpsi_dr) || abs(merluc.dpsi_dr) <= eps
        error('test_salpha_NT:BadFluxDerivative', ...
            'The signed Miller dpsi/dr must be finite and nonzero.');
    end
    Vprime = dVdr/merluc.dpsi_dr;

    if any(~isfinite([V, Vprime, dVdr])) ...
            || V <= 0 || dVdr <= 0 || Vprime == 0
        error('test_salpha_NT:BadVolumeIntegral', ...
            'The Miller V, Vprime, and dV/dr values are invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
