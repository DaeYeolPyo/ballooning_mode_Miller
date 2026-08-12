clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));
addpath(fullfile(this_dir, '..', 'local_equilibrium'));

%% Load equilibrium
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_PT0.6');
%geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_circular');

target_psiN = 0.7;
ntheta0 = 15;
nshat = 50;
shat_end = 5.0;
nalpha = 50;
alpha_end = 6.0;
theta_bnd = 5*pi;

% Use the fine grid for fitting the equilibrium surface, but not for the
% repeated eigenvalue solves.  Using all 1024 surface points on a 10*pi
% ballooning domain creates a dense eigenproblem with about 10,000 DOFs at
% every (s_hat, alpha, theta0) point.
n_surface_points = 1024;
n_balloon_points = 151;

setup_tic = tic;
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);

surf = extract_flux_surface(eq, eqfunc, target_psiN, ...
    Npoints=n_surface_points);
ints = compute_contour_integrals(eq, eqfunc, surf, ...
    isNormalized=false);
[param, bnd] = fit_Miller(eq, eqfunc, target_psiN, ...
    NTheta=n_surface_points);
merluc = Miller_Mercier_Luc(surf, param, bnd);

FFprime = eqfunc.FFprime(target_psiN);
theta0 = linspace(-pi, pi, ntheta0 + 1);
theta0(end) = [];  % -pi and pi represent the same ballooning phase
alpha = linspace(0, alpha_end, nalpha);
shat = linspace(0, shat_end, nshat);
theta_grid = linspace(-theta_bnd, theta_bnd, n_balloon_points).';

aN = eq.rmaxis;
BN = abs(eq.bcentr);

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
    qprime = shat(ss)*(surf.q*ints.Vprime/2/ints.V);

    for aa = 1:nalpha
        pprime = -alpha(aa)*(4*pi^2/2/ints.Vprime)* ...
            sqrt(2*pi^2*eq.rmaxis/ints.V)/(4*pi*1.e-7);

        [~, ~, theta_diag] = evaluate_theta_psi(merluc, ...
            FFprime, pprime, qprime, EnforcePeriodicity=true);

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

figure('Color', 'w', 'Name', 's-alpha stability diagram');
contourf(alpha, shat, lambda_max, 30, 'LineStyle', 'none');
hold on;
contour(alpha, shat, lambda_max, [0, 0], ...
    'k', 'LineWidth', 2);
colorbar;
grid on;
xlabel('$\alpha$', 'Interpreter', 'latex');
ylabel('$\hat{s}$', 'Interpreter', 'latex');

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
cb.Label.String = '$\theta_{0,\mathrm{max}}/\pi$';
cb.Label.Interpreter = 'latex';

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