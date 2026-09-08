clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Indented-Miller equilibrium and scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_I0.8');
target_psiN = 0.7;

ntheta0 = 16;
nshat = 50;
shat_end = 7.0;
nalpha = 50;
alpha_end = 12.0;
theta_bnd = 5*pi;
n_surface_points = int32(1024);
n_fit_points = int32(256);
n_balloon_points = 151;

%% Fit the indented-Miller local equilibrium once
setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=n_surface_points);
ints = compute_contour_integrals(eq, eqfunc, surf, ...
    isNormalized=false);

[param, bnd] = fit_Indented_Miller( ...
    eq, eqfunc, target_psiN, ...
    NTheta=n_fit_points, NRadialFit=int32(9), dpsi=3.e-2, ...
    MaxIterations=int32(250), NPhaseStarts=int32(32), ...
    NXStarts=int32(33));
merluc = Indented_Miller_Mercier_Luc(surf, param, bnd);

[volume, Vprime, dVda] = indented_volume_integrals(merluc);
normal_projection_a = ...
    merluc.normal_sign*merluc.normal_projection_raw;
normal_projection_psi = sign(merluc.dpsi_da)*normal_projection_a;
assert(all(normal_projection_psi > 0), ...
    'The indented-Miller normal is not aligned with grad(psi).');

eq_on_shape = eqfunc.eval(merluc.R, merluc.Z);
Bp_fit_relative_rms = relative_rms(merluc.Bp, eq_on_shape.Bp);
psiN_fit_rms = sqrt(mean((eq_on_shape.psiN(:) - target_psiN).^2));
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

% Definitions shared by the other local-equilibrium s-alpha tests:
%   shat = 2*V*qprime/(q*Vprime)
%   alpha = -mu0*pprime*(2*Vprime)/(4*pi^2)
%           *sqrt(V/(2*pi^2*R_axis)).
qprime_scan = shat.*(surf.q*Vprime/(2*volume));
pprime_scan = -alpha.*(4*pi^2/(2*Vprime)) ...
    .*sqrt(2*pi^2*eq.rmaxis/volume)/mu0;

qprime_equilibrium = eqfunc.qprime(target_psiN);
pprime_equilibrium = eqfunc.pprime(target_psiN);
shat_equilibrium = 2*volume*qprime_equilibrium/(surf.q*Vprime);
alpha_equilibrium = -mu0*pprime_equilibrium*(2*Vprime)/(4*pi^2) ...
    *sqrt(volume/(2*pi^2*eq.rmaxis));

fprintf('Indented-Miller local-equilibrium s-alpha setup\n');
fprintf('  GEQDSK                      : %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  R0, a, delta                : %.8e, %.8e, %+.8e\n', ...
    param.R0, param.a, param.delta);
fprintf('  kappa, I0                   : %.8e, %.8e\n', ...
    param.kappa, param.I0);
fprintf('  dR0/da, dx/da               : %+.8e, %+.8e 1/m\n', ...
    param.dR0_da, param.dx_da);
fprintf('  db/da, dd/da                : %+.8e, %+.8e\n', ...
    param.db_da, param.dd_da);
fprintf('  fitted-path psiN RMS        : %.6e\n', psiN_fit_rms);
fprintf('  Bp relative RMS             : %.6e\n', ...
    Bp_fit_relative_rms);
fprintf('  q target / check            : %.10f / %.10f\n', ...
    surf.q, merluc.q_check);
fprintf('  V fit / GEQDSK              : %.8e / %.8e m^3\n', ...
    volume, ints.V);
fprintf('  |Vprime| fit / GEQDSK       : %.8e / %.8e m^3/Wb\n', ...
    abs(Vprime), abs(ints.Vprime));
fprintf('  V / Vprime relative errors  : %.6e / %.6e\n', ...
    volume_relative_error, Vprime_relative_error);
fprintf('  dV/da, signed dpsi/da       : %.8e, %+.8e\n', ...
    dVda, merluc.dpsi_da);
fprintf('  equilibrium shat / alpha    : %.8f / %.8f\n', ...
    shat_equilibrium, alpha_equilibrium);

assert(Bp_fit_relative_rms < 2.e-2, ...
    'The indented-Miller Bp fit exceeds 2%% RMS error.');
assert(abs(merluc.q_check - surf.q) ...
        < 1.e-10*max(1, abs(surf.q)), ...
    'The indented-Miller local equilibrium does not reproduce q.');
assert(volume_relative_error < 1.e-2 ...
        && Vprime_relative_error < 1.e-2, ...
    'The volume normalization differs from GEQDSK by more than 1%%.');

%% shat -> alpha -> theta0 scan
lambda_max = nan(nshat, nalpha);
best_theta0 = nan(nshat, nalpha);
lambda_theta0 = nan(nshat, nalpha, ntheta0);
fail_message = strings(nshat, nalpha, ntheta0);

nsolve = nshat*nalpha*ntheta0;
fprintf('Equilibrium setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting indented-Miller s-alpha scan: %d shat x ', ...
    '%d alpha x %d theta0 = %d solves.\n'], ...
    nshat, nalpha, ntheta0, nsolve);

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

                % Preserve the secular shear term while reducing the
                % repeated finite-element eigenproblem to a practical grid.
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
fprintf(['Indented-Miller s-alpha scan finished in %.1f min ', ...
    '(%d failed solves).\n'], elapsed/60, nfailed);

%% Package and validate scan products
salpha_Indented_Miller_result = struct();
salpha_Indented_Miller_result.model = ...
    'indented-miller-mercier-luc';
salpha_Indented_Miller_result.shat = shat;
salpha_Indented_Miller_result.alpha = alpha;
salpha_Indented_Miller_result.theta0 = theta0;
salpha_Indented_Miller_result.lambda_max = lambda_max;
salpha_Indented_Miller_result.lambda_theta0 = lambda_theta0;
salpha_Indented_Miller_result.best_theta0 = best_theta0;
salpha_Indented_Miller_result.fail_message = fail_message;
salpha_Indented_Miller_result.qprime = qprime_scan;
salpha_Indented_Miller_result.pprime = pprime_scan;
salpha_Indented_Miller_result.target_psiN = target_psiN;
salpha_Indented_Miller_result.geqdsk_file = geqdsk_file;
salpha_Indented_Miller_result.param = param;
salpha_Indented_Miller_result.volume = volume;
salpha_Indented_Miller_result.Vprime = Vprime;
salpha_Indented_Miller_result.dVda = dVda;
salpha_Indented_Miller_result.volume_relative_error = ...
    volume_relative_error;
salpha_Indented_Miller_result.Vprime_relative_error = ...
    Vprime_relative_error;
salpha_Indented_Miller_result.Bp_fit_relative_rms = ...
    Bp_fit_relative_rms;
salpha_Indented_Miller_result.psiN_fit_rms = psiN_fit_rms;
salpha_Indented_Miller_result.shat_equilibrium = shat_equilibrium;
salpha_Indented_Miller_result.alpha_equilibrium = alpha_equilibrium;
salpha_Indented_Miller_result.scan_elapsed_seconds = elapsed;
salpha_Indented_Miller_result.nfailed = nfailed;

finite_lambda = lambda_max(isfinite(lambda_max));
if isempty(finite_lambda)
    error('test_salpha_Indented_Miller:NoFiniteEigenvalues', ...
        'The s-alpha scan produced no finite eigenvalues.');
end
has_marginal_crossing = min(finite_lambda) <= 0 ...
    && max(finite_lambda) >= 0;
salpha_Indented_Miller_result.has_marginal_crossing = ...
    has_marginal_crossing;

fprintf('  lambda_max range            : [%+.8e, %+.8e]\n', ...
    min(finite_lambda), max(finite_lambda));
fprintf('  marginal crossing found     : %d\n', ...
    has_marginal_crossing);

assert(nfailed == 0, ...
    'The indented-Miller s-alpha scan contains failed solves.');
assert(all(isfinite(lambda_max), 'all'), ...
    'The indented-Miller s-alpha map contains non-finite eigenvalues.');
assert(all(lambda_max(:, 1) < 0), ...
    'The alpha=0 column should be stable on the Dirichlet domain.');

%% Stability diagram and marginal contour
stability_figure = figure('Color', 'w', 'Name', ...
    'Indented-Miller s-alpha stability diagram');
ax_stability = axes(stability_figure);
contourf(ax_stability, alpha, shat, lambda_max, 30, ...
    'LineStyle', 'none');
hold(ax_stability, 'on');
if has_marginal_crossing
    contour(ax_stability, alpha, shat, lambda_max, [0, 0], ...
        'k', 'LineWidth', 2.2);
else
    warning('test_salpha_Indented_Miller:NoMarginalCrossing', ...
        ['No lambda_max=0 crossing was found. Increase shat_end or ', ...
         'alpha_end, or refine the scan range.']);
end

equilibrium_is_visible = alpha_equilibrium >= min(alpha) ...
    && alpha_equilibrium <= max(alpha) ...
    && shat_equilibrium >= min(shat) ...
    && shat_equilibrium <= max(shat);
if ~equilibrium_is_visible
    plot(ax_stability, alpha_equilibrium, shat_equilibrium, 'p', ...
        'MarkerSize', 13, 'MarkerEdgeColor', 'w', ...
        'MarkerFaceColor', [0.85, 0.1, 0.75], 'LineWidth', 1.4);
end

cb_stability = colorbar(ax_stability);
cb_stability.Label.Interpreter = 'latex';
cb_stability.Label.String = '$\lambda_{\max}$';
grid(ax_stability, 'on');
box(ax_stability, 'on');
ax_stability.Layer = 'top';
ax_stability.TickDir = 'out';
xlabel(ax_stability, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax_stability, '$\hat{s}$', 'Interpreter', 'latex');

%% Ballooning phase that maximizes lambda
[alpha_grid, shat_grid] = meshgrid(alpha, shat);
valid = isfinite(best_theta0);

phase_figure = figure('Color', 'w', 'Name', ...
    'Indented-Miller phase of maximum eigenvalue');
ax_phase = axes(phase_figure);
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
cb_phase.Label.String = '$\theta_{0,\max}/\pi$';

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
title(ax_phase, ...
    'Indented-Miller phase of maximum ballooning eigenvalue');

salpha_Indented_Miller_result.stability_figure = stability_figure;
salpha_Indented_Miller_result.phase_figure = phase_figure;
salpha_Indented_Miller_result.status = 'PASS';

fprintf('  validation status           : PASS\n');

function [V, Vprime, dVda] = indented_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVda = 2*pi*trapz(theta, R.*abs(jac));

    if ~isfinite(merluc.dpsi_da) || abs(merluc.dpsi_da) <= eps
        error('test_salpha_Indented_Miller:BadFluxDerivative', ...
            'The signed indented-Miller dpsi/da must be nonzero.');
    end
    Vprime = dVda/merluc.dpsi_da;

    if any(~isfinite([V, Vprime, dVda])) ...
            || V <= 0 || dVda <= 0 || Vprime == 0
        error('test_salpha_Indented_Miller:BadVolumeIntegral', ...
            'The fitted V, Vprime, or dV/da value is invalid.');
    end
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
