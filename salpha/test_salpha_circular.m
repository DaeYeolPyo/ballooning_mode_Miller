clc
close all
clear

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Circular equilibrium and scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_circular');

if ~exist('target_psiN', 'var') || isempty(target_psiN)
    target_psiN = 0.8;
end
if ~exist('ntheta0', 'var') || isempty(ntheta0)
    ntheta0 = 16;
end
if ~exist('nshat', 'var') || isempty(nshat)
    nshat = 50;
end
if ~exist('shat_end', 'var') || isempty(shat_end)
    shat_end = 5.0;
end
if ~exist('nalpha', 'var') || isempty(nalpha)
    nalpha = 50;
end
if ~exist('alpha_end', 'var') || isempty(alpha_end)
    alpha_end = 6.0;
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

%% Fit the circular local equilibrium once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=int32(n_surface_points));

[param, bnd] = fit_Circular(eq, eqfunc, target_psiN, ...
    NTheta=int32(n_surface_points));
circluc = Circular_Mercier_Luc(surf, param, bnd);

[volume, Vprime, dVdr, volume_numeric_error] = ...
    circular_volume_integrals(circluc, param);

normal_projection_r = ...
    circluc.normal_sign*circluc.normal_projection_raw;
normal_projection_psi = sign(circluc.dpsi_dr)*normal_projection_r;
assert(all(normal_projection_psi > 0), ...
    'The circular Mercier normal is not aligned with grad(psi).');

eq_on_circle = eqfunc.eval(circluc.R, circluc.Z);
Bp_fit_relative_rms = relative_rms(circluc.Bp, eq_on_circle.Bp);
psiN_fit_rms = sqrt(mean((eq_on_circle.psiN(:) - target_psiN).^2));

FFprime = eqfunc.FFprime(target_psiN);
% Uniform cyclic phases with theta0=0 included for both odd and even counts.
theta0 = 2*pi*((0:ntheta0-1) - floor(ntheta0/2))/ntheta0;
alpha = linspace(0, alpha_end, nalpha);
shat = linspace(0, shat_end, nshat);
theta_grid = linspace(-theta_bnd, theta_bnd, n_balloon_points).';

aN = eq.rmaxis;
BN = abs(eq.bcentr);
mu0 = 4*pi*1.e-7;

fprintf('Circular local-equilibrium s-alpha setup\n');
fprintf('  GEQDSK                      : %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  R0, r, dR0/dr              : %.8e, %.8e, %.8e\n', ...
    param.R0, param.r, param.dR0_dr);
fprintf('  contour radial RMS / max   : %.6e / %.6e m\n', ...
    param.fit_residual_rms, param.fit_residual_max);
fprintf('  circular-path psiN RMS     : %.6e\n', psiN_fit_rms);
fprintf('  Bp relative RMS            : %.6e\n', Bp_fit_relative_rms);
fprintf('  q target / check           : %.10f / %.10f\n', ...
    surf.q, circluc.q_check);
fprintf('  analytic/numeric q-int err : %.3e\n', ...
    circluc.q_integral_relative_error);
fprintf('  V, dV/dr, Vprime           : %.8e, %.8e, %+.8e\n', ...
    volume, dVdr, Vprime);
fprintf('  analytic/numeric volume err: %.3e\n', volume_numeric_error);

assert(Bp_fit_relative_rms < 5.e-2, ...
    'The circular Bp fit exceeds the 5%% RMS validation threshold.');
assert(abs(circluc.q_check - surf.q) < 1.e-10*max(1, abs(surf.q)), ...
    'The circular local equilibrium does not reproduce q.');
assert(circluc.q_integral_relative_error < 1.e-10, ...
    'The circular analytic q integral is inconsistent.');
assert(volume_numeric_error < 1.e-10, ...
    'The circular analytic and numerical volumes are inconsistent.');

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
fprintf(['Starting circular s-alpha scan: %d shat x %d alpha x ', ...
    '%d theta0 = %d solves.\n'], nshat, nalpha, ntheta0, nsolve);

scan_tic = tic;
for ss = 1:nshat
    row_tic = tic;

    % Definitions used by the other local-equilibrium s-alpha tests:
    %   shat = 2*V*qprime/(q*Vprime)
    %   alpha = -mu0*pprime*(2*Vprime)/(4*pi^2)
    %           *sqrt(V/(2*pi^2*R_axis)).
    qprime = qprime_scan(ss);

    for aa = 1:nalpha
        pprime = pprime_scan(aa);

        [~, ~, theta_diag] = evaluate_theta_psi(circluc, ...
            FFprime, pprime, qprime, EnforcePeriodicity=true);

        for tt = 1:ntheta0
            try
                bal = evaluate_gcf(theta_diag, theta0(tt), ...
                    theta_bnd=theta_bnd, a_N=aN, B_N=BN);

                % Preserve the non-periodic shear term while reducing the
                % finite-element problem to a practical scan grid.
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
fprintf('Circular s-alpha scan finished in %.1f min (%d failed solves).\n', ...
    elapsed/60, nfailed);

%% Package scan products
salpha_circular_result = struct();
salpha_circular_result.model = 'circular-mercier-luc';
salpha_circular_result.shat = shat;
salpha_circular_result.alpha = alpha;
salpha_circular_result.theta0 = theta0;
salpha_circular_result.lambda_max = lambda_max;
salpha_circular_result.lambda_theta0 = lambda_theta0;
salpha_circular_result.best_theta0 = best_theta0;
salpha_circular_result.fail_message = fail_message;
salpha_circular_result.qprime = qprime_scan;
salpha_circular_result.pprime = pprime_scan;
salpha_circular_result.target_psiN = target_psiN;
salpha_circular_result.geqdsk_file = geqdsk_file;
salpha_circular_result.param = param;
salpha_circular_result.volume = volume;
salpha_circular_result.Vprime = Vprime;
salpha_circular_result.dVdr = dVdr;
salpha_circular_result.Bp_fit_relative_rms = Bp_fit_relative_rms;
salpha_circular_result.psiN_fit_rms = psiN_fit_rms;
salpha_circular_result.scan_elapsed_seconds = elapsed;
salpha_circular_result.nfailed = nfailed;

finite_lambda = lambda_max(isfinite(lambda_max));
if isempty(finite_lambda)
    error('test_salpha_circular:NoFiniteEigenvalues', ...
        'The circular s-alpha scan produced no finite eigenvalues.');
end
has_marginal_crossing = min(finite_lambda) <= 0 ...
    && max(finite_lambda) >= 0;
salpha_circular_result.has_marginal_crossing = has_marginal_crossing;

fprintf('  lambda_max range            : [%+.8e, %+.8e]\n', ...
    min(finite_lambda), max(finite_lambda));
fprintf('  marginal crossing found     : %d\n', has_marginal_crossing);

assert(nfailed == 0, ...
    'The circular s-alpha scan contains failed eigenvalue solves.');
assert(all(isfinite(lambda_max), 'all'), ...
    'The circular s-alpha map contains non-finite eigenvalues.');
assert(all(lambda_max(:, 1) < 0), ...
    'The alpha=0 column should be stable on the finite Dirichlet domain.');

%% Stability and marginal-contour plot
figure('Color', 'w', 'Name', 'Circular s-alpha stability diagram');
ax_stability = axes();
contourf(ax_stability, alpha, shat, lambda_max, 30, ...
    'LineStyle', 'none');
hold(ax_stability, 'on');
if has_marginal_crossing
    contour(ax_stability, alpha, shat, lambda_max, [0, 0], ...
        'k', 'LineWidth', 2.2);
else
    warning('test_salpha_circular:NoMarginalCrossing', ...
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
title(ax_stability, sprintf(['Circular local-equilibrium $s$-$\\alpha$ ', ...
    'stability, $\\psi_N=%.2f$'], target_psiN), ...
    'Interpreter', 'latex');

%% Ballooning phase that maximizes lambda
[alpha_grid, shat_grid] = meshgrid(alpha, shat);
valid = isfinite(best_theta0);

figure('Color', 'w', 'Name', ...
    'Circular phase of maximum ballooning eigenvalue');
ax_phase = axes();
scatter(ax_phase, alpha_grid(valid), shat_grid(valid), 58, ...
    best_theta0(valid)/pi, 's', 'filled', ...
    'MarkerEdgeColor', [0.18, 0.18, 0.18], 'LineWidth', 0.25);

theta0_over_pi = theta0/pi;
dtheta0_over_pi = 2/ntheta0;
colormap(ax_phase, hsv(ntheta0));
clim(ax_phase, [theta0_over_pi(1) - dtheta0_over_pi/2, ...
                theta0_over_pi(end) + dtheta0_over_pi/2]);
cb = colorbar(ax_phase);
cb.Ticks = theta0_over_pi;
cb.TickLabels = compose('%+.2f', theta0_over_pi);
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\theta_{0,\mathrm{max}}/\pi$';

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
title(ax_phase, 'Circular phase of maximum ballooning eigenvalue');

fprintf('  validation status           : PASS\n');

function [V, Vprime, dVdr, numeric_error] = ...
        circular_volume_integrals(circluc, param)
    theta = circluc.theta_geo(:);
    R = circluc.R(:);
    Zt = circluc.Zt(:);
    jac = circluc.jac(:);

    % Shifted circular torus:
    %   V = 2*pi^2*R0*r^2,
    %   dV/dr = 2*pi^2*(2*R0*r + dR0/dr*r^2).
    V = 2*pi^2*param.R0*param.r^2;
    dVdr = 2*pi^2*( ...
        2*param.R0*param.r + param.dR0_dr*param.r^2);

    V_numeric = pi*abs(trapz(theta, R.^2.*Zt));
    dVdr_numeric = 2*pi*trapz(theta, R.*abs(jac));
    numeric_error = max([ ...
        abs(V_numeric - V)/max(abs(V), eps), ...
        abs(dVdr_numeric - dVdr)/max(abs(dVdr), eps)]);

    if ~isfinite(circluc.dpsi_dr) || abs(circluc.dpsi_dr) <= eps
        error('test_salpha_circular:BadFluxDerivative', ...
            'The circular dpsi/dr must be finite and nonzero.');
    end
    Vprime = dVdr/circluc.dpsi_dr;

    if any(~isfinite([V, Vprime, dVdr])) ...
            || V <= 0 || dVdr <= 0 || Vprime == 0
        error('test_salpha_circular:BadVolumeIntegral', ...
            'The circular V, Vprime, and dV/dr values are invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
