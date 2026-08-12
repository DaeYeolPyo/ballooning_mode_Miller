clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Load equilibrium and set scan controls
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'kappa5', ...
    ['scan_B2.5_C10_G13_H12_app2.033_gpp0.935_', ...
     'aff0.938_gff1.800.geqdsk']);

target_psiN = 0.7;
ntheta0 = 15;
nshat = 50;
shat_end = 5.0;
nalpha = 50;
alpha_end = 6.0;
theta_bnd = 5*pi;

% The C-shape fit is converged primarily by its alternating iterations; 256
% geometric points are sufficient.  The repeated eigenvalue solves use a
% separate, deliberately smaller ballooning grid.
n_surface_points = int32(1024);
n_cshape_points = int32(256);
n_radial_fit = int32(9);
max_fit_iterations = int32(3000);
radial_fit_half_width = 3.e-2;
n_balloon_points = 151;

setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);

surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=n_surface_points);
ints = compute_contour_integrals(eq, eqfunc, surf, ...
    isNormalized=false);
[param, bnd] = fit_CShape(eq, eqfunc, target_psiN, ...
    NTheta=n_cshape_points, ...
    NRadialFit=n_radial_fit, ...
    dpsi=radial_fit_half_width, ...
    MaxIterations=max_fit_iterations);
merluc = CShape_Mercier_Luc(surf, param, bnd);

% Use the local C-shape volume derivatives when converting the requested
% s-alpha coordinates to qprime and pprime.  Vprime remains signed because
% evaluate_theta_psi differentiates with respect to physical GEQDSK psi.
[volume, Vprime, dVdr] = cshape_volume_integrals(merluc);

normal_projection_r = ...
    merluc.normal_sign*merluc.normal_projection_raw;
normal_projection_psi = sign(merluc.dpsi_dr)*normal_projection_r;
assert(all(normal_projection_psi > 0), ...
    'The C-shape Mercier normal is not aligned with grad(psi).');

% CShape_boundary starts near the upper point.  Report theta0 relative to
% the conventional outboard PEST origin used by the direct equilibrium.
[~, outboard_index] = max(merluc.R(1:end-1));
theta_origin_shift = merluc.theta_PEST(outboard_index);

FFprime = eqfunc.FFprime(target_psiN);
theta0 = linspace(-pi, pi, ntheta0 + 1);
theta0(end) = [];  % -pi and pi represent the same ballooning phase
alpha = linspace(0, alpha_end, nalpha);
shat = linspace(0, shat_end, nshat);
theta_grid = linspace(-theta_bnd, theta_bnd, n_balloon_points).';

aN = eq.rmaxis;
BN = abs(eq.bcentr);
mu0 = 4*pi*1.e-7;

eq_on_cshape = eqfunc.eval(merluc.R, merluc.Z);
Bp_fit_relative_rms = relative_rms(merluc.Bp, eq_on_cshape.Bp);

fprintf('C-shape local-equilibrium setup\n');
fprintf('  target psiN                : %.6f\n', target_psiN);
fprintf('  coefficients [A B C G H]  :');
fprintf(' %+.8e', param.coeff);
fprintf('\n');
fprintf('  radial derivatives         :');
fprintf(' %+.8e', param.dcoeff_dr);
fprintf('\n');
fprintf('  fit relative RMS           : %.6e\n', ...
    param.fit.relative_rms_distance);
fprintf('  Bp relative RMS            : %.6e\n', Bp_fit_relative_rms);
fprintf('  q target / check           : %.10f / %.10f\n', ...
    surf.q, merluc.q_check);
fprintf('  V C-shape / GEQDSK         : %.8e / %.8e m^3\n', ...
    volume, ints.V);
fprintf('  |Vprime| C-shape / GEQDSK  : %.8e / %.8e m^3/Wb\n', ...
    abs(Vprime), ints.Vprime);
fprintf('  dV/dr, signed dpsi/dr      : %.8e, %+.8e\n', ...
    dVdr, merluc.dpsi_dr);
fprintf('  outboard theta-origin shift: %.8e rad\n', ...
    theta_origin_shift);

lambda_max = nan(nshat, nalpha);
best_theta0 = nan(nshat, nalpha);
lambda_theta0 = nan(nshat, nalpha, ntheta0);
fail_message = strings(nshat, nalpha, ntheta0);

nsolve = nshat*nalpha*ntheta0;
fprintf('Equilibrium setup finished in %.2f s.\n', toc(setup_tic));
fprintf(['Starting s-alpha scan: %d shat x %d alpha x %d theta0 ', ...
    '= %d solves.\n'], nshat, nalpha, ntheta0, nsolve);

scan_tic = tic;
for ss = 1:nshat
    row_tic = tic;
    qprime = shat(ss)*(surf.q*Vprime/(2*volume));

    for aa = 1:nalpha
        pprime = -alpha(aa)*(4*pi^2/(2*Vprime))* ...
            sqrt(2*pi^2*eq.rmaxis/volume)/mu0;

        [~, ~, theta_diag] = evaluate_theta_psi(merluc, ...
            FFprime, pprime, qprime, EnforcePeriodicity=true);
        theta_diag = shift_theta_origin(theta_diag, theta_origin_shift);

        for tt = 1:ntheta0
            t0 = theta0(tt);

            try
                bal = evaluate_gcf(theta_diag, t0, ...
                    theta_bnd=theta_bnd, a_N=aN, B_N=BN);

                % Resample the already-extended, non-periodic coefficients.
                % This preserves the secular shear term while keeping the
                % finite-element eigenproblem small enough for a scan.
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
    fprintf(['  shat row %2d/%2d done: shat = %.4g, ', ...
        'row %.1f s, ETA %.1f min\n'], ...
        ss, nshat, shat(ss), toc(row_tic), eta/60);
end

elapsed = toc(scan_tic);
nfailed = nnz(strlength(fail_message) > 0);
fprintf('s-alpha scan finished in %.1f min (%d failed solves).\n', ...
    elapsed/60, nfailed);

% Keep the scan products together for later plotting or saving.  Positive
% lambda_max is unstable; retaining negative values is necessary to draw
% the marginal-stability contour lambda_max = 0.
salpha_result = struct();
salpha_result.shat = shat;
salpha_result.alpha = alpha;
salpha_result.theta0 = theta0;
salpha_result.lambda_max = lambda_max;
salpha_result.lambda_theta0 = lambda_theta0;
salpha_result.best_theta0 = best_theta0;
salpha_result.fail_message = fail_message;
salpha_result.target_psiN = target_psiN;
salpha_result.geqdsk_file = geqdsk_file;
salpha_result.model = 'cshape-mercier-luc';
salpha_result.param = param;
salpha_result.volume = volume;
salpha_result.Vprime = Vprime;
salpha_result.dVdr = dVdr;
salpha_result.theta_origin_shift = theta_origin_shift;
salpha_result.Bp_fit_relative_rms = Bp_fit_relative_rms;
%%

figure('Color', 'w', 'Name', 's-alpha stability diagram');
contourf(alpha, shat, lambda_max, 30, 'LineStyle', 'none');
hold on;
contour(alpha, shat, lambda_max, [0, 0], ...
    'k', 'LineWidth', 2);
colorbar;
grid on;
xlabel('$\alpha$', 'Interpreter', 'latex');
ylabel('$\hat{s}$', 'Interpreter', 'latex');
title('C-shape fitting based $s$-$\alpha$ stability', ...
    'Interpreter', 'latex');

% best_theta0 only takes values from the discrete theta0 scan.  Plot one
% colored square per (s_hat, alpha) sample instead of interpolating between
% grid points with contourf.
[alpha_grid, shat_grid] = meshgrid(alpha, shat);
valid = isfinite(best_theta0);

figure('Color', 'w', 'Name', 'Phase of maximum ballooning eigenvalue');
ax = axes();
scatter(ax, alpha_grid(valid), shat_grid(valid), 58, ...
    best_theta0(valid)/pi, 's', 'filled', ...
    'MarkerEdgeColor', [0.18, 0.18, 0.18], 'LineWidth', 0.25);

% A cyclic map is appropriate because theta0 = -pi and theta0 = pi are
% the same ballooning phase.  The half-step color limits make each sampled
% theta0 value occupy exactly one discrete color band.
theta0_over_pi = theta0/pi;
dtheta0_over_pi = 2/ntheta0;
colormap(ax, hsv(ntheta0));
clim(ax, [theta0_over_pi(1) - dtheta0_over_pi/2, ...
          theta0_over_pi(end) + dtheta0_over_pi/2]);

cb = colorbar(ax);
cb.Ticks = theta0_over_pi;
cb.TickLabels = compose('%+.2f', theta0_over_pi);
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\theta_{0,\mathrm{max}}/\pi$';

if numel(alpha) > 1
    dalpha = alpha(2) - alpha(1);
else
    dalpha = 1;
end
if numel(shat) > 1
    dshat = shat(2) - shat(1);
else
    dshat = 1;
end
xlim(ax, [min(alpha) - 0.5*dalpha, max(alpha) + 0.5*dalpha]);
ylim(ax, [min(shat) - 0.5*dshat, max(shat) + 0.5*dshat]);

grid(ax, 'on');
box(ax, 'on');
ax.Layer = 'top';
ax.TickDir = 'out';
ax.FontSize = 11;
xlabel(ax, '$\alpha$', 'Interpreter', 'latex');
ylabel(ax, '$\hat{s}$', 'Interpreter', 'latex');
title(ax, 'C-shape phase of maximum ballooning eigenvalue');

function [V, Vprime, dVdr] = cshape_volume_integrals(merluc)
    theta = merluc.theta_geo(:);
    R = merluc.R(:);
    Zt = merluc.Zt(:);
    jac = merluc.jac(:);

    % V is the enclosed toroidal volume.  r=B is the increasing-outward
    % C-shape radial coordinate, so dV/dr is positive.  Vprime=dV/dpsi is
    % signed and follows the physical GEQDSK flux convention.
    V = pi*abs(trapz(theta, R.^2.*Zt));
    dVdr = 2*pi*trapz(theta, R.*abs(jac));

    if ~isfinite(merluc.dpsi_dr) || abs(merluc.dpsi_dr) <= eps
        error('test_salpha_Cshape:BadFluxDerivative', ...
            'The signed C-shape dpsi/dr must be finite and nonzero.');
    end
    Vprime = dVdr/merluc.dpsi_dr;

    if any(~isfinite([V, Vprime, dVdr])) ...
            || V <= 0 || dVdr <= 0 || Vprime == 0
        error('test_salpha_Cshape:BadVolumeIntegral', ...
            'The C-shape volume derivatives are invalid.');
    end
end

function diag = shift_theta_origin(diag, theta_shift)
    theta = diag.theta(:);
    theta_query = theta + theta_shift;

    periodic_fields = { ...
        'R', 'Z', 'Bp', 'Bphi', 'B2', 'B', 'H', 'dl_dt', ...
        'normal_R', 'normal_Z', 'tangent_R', 'tangent_Z', ...
        'curvature_normal', 'dBp_dl', 'Gamma'};

    for k = 1:numel(periodic_fields)
        name = periodic_fields{k};
        if isfield(diag, name)
            diag.(name) = interp_periodic_closed( ...
                theta, diag.(name), theta_query);
        end
    end

    theta_psi_shifted = interp_periodic_closed( ...
        theta, diag.theta_psi, theta_query);
    diag.theta_psi = theta_psi_shifted - theta_psi_shifted(1);
    diag.theta_origin_shift = theta_shift;
end

function yq = interp_periodic_closed(theta, y, thetaq)
    theta = theta(:);
    y = y(:);
    theta_core = theta(1:end-1);
    y_core = y(1:end-1);
    period = theta(end) - theta(1);
    theta_ext = [theta_core - period; theta_core; theta_core + period];
    y_ext = [y_core; y_core; y_core];
    theta_wrapped = mod(thetaq(:) - theta(1), period) + theta(1);
    yq = interp1(theta_ext, y_ext, theta_wrapped, 'pchip');
end

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end
