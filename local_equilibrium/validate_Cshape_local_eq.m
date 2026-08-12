clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));
addpath(fullfile(this_dir, '..', 'ballooning_equation'));
addpath(fullfile(this_dir, '..', 'sturm_liouville_solver'));

%% Common equilibrium, normalization, and ballooning domain
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'kappa5', ...
    ['scan_B2.5_C10_G13_H12_app2.033_gpp0.935_', ...
     'aff0.938_gff1.800.geqdsk']);

target_psiN = 0.7;
theta0 = 0.0;
theta_bnd = 5*pi;

n_theta_metric = 256;
n_theta_compare = 401;
n_surface_points = int32(1024);

% C-shape fitting controls.  NTheta=256 is already sufficient for the
% boundary representation; MaxIterations=1000 gives reasonably converged
% coefficients without making this validation prohibitively expensive.
n_theta_cshape = int32(256);
n_radial_fit = int32(9);
max_fit_iterations = int32(1000);
radial_fit_half_width = 3.e-2;

eq_physical = read_geqdsk(geqdsk_file);
eqfunc_physical = build_interpolants(eq_physical);

% scale_eq uses these normalizations.  The C-shape path must use the same
% a_N and B_N for its dimensionless g, c, and f to be directly comparable.
a_N = eq_physical.rmaxis;
B_N = abs(eq_physical.bcentr);
psi_scale = a_N^2*B_N;

theta_common = linspace(-theta_bnd, theta_bnd, n_theta_compare).';

%% Path A: direct PEST metric from the Grad-Shafranov equilibrium
eq_metric = scale_eq(eq_physical);
eqfunc_metric = build_interpolants(eq_metric);

radial_step = 5.e-3;
psiN_list = target_psiN + (-4:4).'*radial_step;
surfaces = repmat(struct('surf', [], 'ints', [], 'ang', []), ...
    numel(psiN_list), 1);

% compute_metrics uses the signed physical-flux coordinate and requires a
% positive (psi,theta,zeta) Jacobian.  extract_flux_surface always returns
% the outboard-up (CCW) ordering, so reverse it when physical psi decreases
% outward.  This preserves the signed Jacobian instead of replacing it by
% its absolute value.
reverse_direct_surfaces = eq_metric.sibry - eq_metric.simag < 0;

fprintf('Building direct PEST metric from %d GEQDSK surfaces ...\n', ...
    numel(psiN_list));
for k = 1:numel(psiN_list)
    surf_k = extract_flux_surface( ...
        eq_metric, eqfunc_metric, psiN_list(k), ...
        Npoints=n_surface_points);
    if reverse_direct_surfaces
        surf_k = reverse_closed_surface(surf_k);
    end
    ints_k = compute_contour_integrals( ...
        eq_metric, eqfunc_metric, surf_k, isNormalized=true);
    ang_k = straight_field_line_angles( ...
        eqfunc_metric, surf_k, ints_k, UseQ="geqdsk");

    surfaces(k).surf = surf_k;
    surfaces(k).ints = ints_k;
    surfaces(k).ang = ang_k;
end

maps = build_SFL_maps(surfaces, n_theta_metric);
map = maps.PEST;
metrics = compute_metrics(map, eqfunc_metric);

coeff_metric = calculate_ballooning_coeffs( ...
    map, metrics, target_psiN, ...
    Theta0=theta0, ThetaExt=theta_common, use_salpha=false);

[~, target_index] = min(abs(map.psiN - target_psiN));
theta_base = map.theta(:);
R_metric = map.R(target_index, :).';
Z_metric = map.Z(target_index, :).';
F_metric = map.F(target_index);
B2_metric = metrics.B2(target_index, :).';
Bp_metric = sqrt(max(B2_metric - (F_metric./R_metric).^2, 0));
theta_psi_metric = metrics.g_contra.psitheta(target_index, :).' ...
    ./metrics.g_contra.psipsi(target_index, :).';

%% Path B: fitted C-shape with Mercier-Luc local equilibrium
surf_local = extract_flux_surface( ...
    eq_physical, eqfunc_physical, target_psiN, ...
    Npoints=n_surface_points);

fprintf('Fitting C-shape on %d radial surfaces ...\n', n_radial_fit);
[param, bnd] = fit_CShape( ...
    eq_physical, eqfunc_physical, target_psiN, ...
    NTheta=n_theta_cshape, ...
    NRadialFit=n_radial_fit, ...
    dpsi=radial_fit_half_width, ...
    MaxIterations=max_fit_iterations);

merluc = CShape_Mercier_Luc(surf_local, param, bnd);

FFprime = eqfunc_physical.FFprime(target_psiN);
pprime = eqfunc_physical.pprime(target_psiN);
dpsiN_q = 0.01;
dq_dpsiN = (eqfunc_physical.q(target_psiN + dpsiN_q) ...
    - eqfunc_physical.q(target_psiN - dpsiN_q))/(2*dpsiN_q);
qprime = dq_dpsiN/(eq_physical.sibry - eq_physical.simag);

[~, ~, theta_diag] = evaluate_theta_psi( ...
    merluc, FFprime, pprime, qprime, EnforcePeriodicity=true);

% The analytic C-shape starts at its upper point, whereas the direct map
% defines theta=0 at the outboard midplane on every surface.  A poloidal
% angle has a flux-function gauge freedom, so rotate all periodic local
% fields to the same outboard origin and set theta_psi=0 there before making
% pointwise or fixed-theta0 comparisons.
[~, outboard_index] = max(merluc.R(1:end-1));
theta_origin_shift = theta_diag.theta(outboard_index);
theta_diag = shift_theta_origin(theta_diag, theta_origin_shift);

bal_local = evaluate_gcf(theta_diag, theta0, ...
    theta_bnd=theta_bnd, a_N=a_N, B_N=B_N);

% The explicit shear dependence is not periodic on the extended ballooning
% domain, so use ordinary rather than periodic interpolation here.
g_local = interp1(bal_local.theta, bal_local.g, theta_common, 'pchip');
c_local = interp1(bal_local.theta, bal_local.c, theta_common, 'pchip');
f_local = interp1(bal_local.theta, bal_local.f, theta_common, 'pchip');
K_local = interp1(bal_local.theta, ...
    a_N^2.*bal_local.K, theta_common, 'pchip');

% Compare one-period geometry on the direct metric's uniform PEST grid.
R_local = interp_periodic_closed( ...
    theta_diag.theta, theta_diag.R./a_N, theta_base);
Z_local = interp_periodic_closed( ...
    theta_diag.theta, theta_diag.Z./a_N, theta_base);
Bp_local = interp_periodic_closed( ...
    theta_diag.theta, theta_diag.Bp./B_N, theta_base);
theta_psi_local_N = interp_periodic_closed( ...
    theta_diag.theta, psi_scale.*theta_diag.theta_psi, theta_base);

% This comparison isolates the C-shape Bp reconstruction on its own fitted
% boundary, independently of the PEST remapping used above.
eq_on_cshape = eqfunc_physical.eval(merluc.R, merluc.Z);
Bp_fit_reference = eq_on_cshape.Bp(:);
Bp_fit_relative_rms = relative_rms(merluc.Bp, Bp_fit_reference);

%% Common finite-element ballooning solve
[G_metric, C_metric, F_metric_mat, theta_dof] = construct_matrix( ...
    theta_common, coeff_metric.g, coeff_metric.c, coeff_metric.f);
[lambda_metric_all, X_metric_all] = calculate_eigenvalue( ...
    theta_dof, G_metric, C_metric, F_metric_mat, 'dirichlet');

[G_local, C_local, F_local_mat, theta_dof_local] = construct_matrix( ...
    theta_common, g_local, c_local, f_local);
[lambda_local_all, X_local_all] = calculate_eigenvalue( ...
    theta_dof_local, G_local, C_local, F_local_mat, 'dirichlet');

lambda_metric = lambda_metric_all(1);
lambda_local = lambda_local_all(1);
mode_metric = normalize_mode(X_metric_all(:, 1));
mode_local = normalize_mode(X_local_all(:, 1));
if trapz(theta_dof, mode_metric.*mode_local) < 0
    mode_local = -mode_local;
end
mode_overlap = abs(trapz(theta_dof, mode_metric.*mode_local)) ...
    /sqrt(trapz(theta_dof, mode_metric.^2) ...
          *trapz(theta_dof, mode_local.^2));

%% Quantitative comparison
comparison = struct();

comparison.geometry.R_relative_rms = relative_rms(R_local, R_metric);
comparison.geometry.Z_absolute_rms = absolute_rms(Z_local, Z_metric);
comparison.geometry.Bp_relative_rms = relative_rms(Bp_local, Bp_metric);
comparison.geometry.Bp_fit_relative_rms = Bp_fit_relative_rms;
comparison.geometry.theta_psi_relative_rms = relative_rms( ...
    theta_psi_local_N, theta_psi_metric);

comparison.profile.q_metric = coeff_metric.q;
comparison.profile.q_local = theta_diag.q;
comparison.profile.qprime_metric = coeff_metric.qprime;
comparison.profile.qprime_local = psi_scale*qprime;
comparison.profile.pprime_metric = coeff_metric.pprime;
comparison.profile.pprime_local = bal_local.dpbar_dpsiN;

comparison.local.FFprime_input = theta_diag.FFprime_input;
comparison.local.FFprime_periodic = theta_diag.FFprime_periodic;
comparison.local.FFprime_relative_difference = abs( ...
    theta_diag.FFprime_periodic - theta_diag.FFprime_input) ...
    /max(abs(theta_diag.FFprime_input), eps);
comparison.local.raw_closure_defect = theta_diag.input_closure_defect;
comparison.local.dpsi_dr_q = merluc.dpsi_dr;
comparison.local.dpsi_dr_fit = param.dpsi_dr_fit;
comparison.local.dpsi_dr_relative_difference = abs( ...
    merluc.dpsi_dr - param.dpsi_dr_fit) ...
    /max(abs(param.dpsi_dr_fit), eps);
comparison.local.theta_PEST_span = ...
    merluc.theta_PEST(end) - merluc.theta_PEST(1);
comparison.local.theta_origin_shift = theta_origin_shift;

normal_projection_r = ...
    merluc.normal_sign*merluc.normal_projection_raw;
normal_projection_psi = sign(merluc.dpsi_dr)*normal_projection_r;
comparison.local.minimum_flux_normal_projection = ...
    min(normal_projection_psi);

comparison.coefficient.g_relative_rms = relative_rms( ...
    g_local, coeff_metric.g);
comparison.coefficient.c_relative_rms = relative_rms( ...
    c_local, coeff_metric.c);
comparison.coefficient.f_relative_rms = relative_rms( ...
    f_local, coeff_metric.f);
comparison.coefficient.K_relative_rms = relative_rms( ...
    K_local, coeff_metric.K);

comparison.eigenvalue.metric = lambda_metric;
comparison.eigenvalue.local = lambda_local;
comparison.eigenvalue.absolute_difference = abs(lambda_local - lambda_metric);
comparison.eigenvalue.relative_difference = abs(lambda_local - lambda_metric) ...
    /max(abs(lambda_metric), eps);
comparison.eigenfunction.overlap = mode_overlap;

fprintf('\nDirect GS metric vs C-shape Mercier-Luc equilibrium\n');
fprintf('  GEQDSK / target psiN       : %s / %.6f\n', ...
    geqdsk_file, target_psiN);
fprintf('  a_N, B_N                   : %.8e m, %.8e T\n', a_N, B_N);
fprintf('  q direct/local             : %.10f / %.10f\n', ...
    comparison.profile.q_metric, comparison.profile.q_local);
fprintf('  qprime_N direct/local      : %.8e / %.8e\n', ...
    comparison.profile.qprime_metric, comparison.profile.qprime_local);
fprintf('  pprime_N direct/local      : %.8e / %.8e\n', ...
    comparison.profile.pprime_metric, comparison.profile.pprime_local);

fprintf('\nC-shape local-equilibrium consistency\n');
fprintf('  dpsi/dr q / radial fit     : %+.8e / %+.8e Wb/m\n', ...
    comparison.local.dpsi_dr_q, comparison.local.dpsi_dr_fit);
fprintf('  relative dpsi/dr difference: %.6e\n', ...
    comparison.local.dpsi_dr_relative_difference);
fprintf('  min flux-normal projection : %.8e\n', ...
    comparison.local.minimum_flux_normal_projection);
fprintf('  theta_PEST span / 2pi      : %.12f / %.12f\n', ...
    comparison.local.theta_PEST_span, 2*pi);
fprintf('  outboard theta-origin shift: %.8e rad\n', ...
    comparison.local.theta_origin_shift);
fprintf('  FFprime input / periodic   : %+.8e / %+.8e\n', ...
    comparison.local.FFprime_input, ...
    comparison.local.FFprime_periodic);
fprintf('  relative FFprime difference: %.6e\n', ...
    comparison.local.FFprime_relative_difference);
fprintf('  raw theta_psi closure defect: %.6e\n', ...
    comparison.local.raw_closure_defect);
fprintf('  vector identity K/D/g/c/f  : %.3e %.3e %.3e %.3e %.3e\n', ...
    bal_local.validation.K_relative_error, ...
    bal_local.validation.D_relative_error, ...
    bal_local.validation.g_relative_error, ...
    bal_local.validation.c_relative_error, ...
    bal_local.validation.f_relative_error);

fprintf('\nOne-period geometry comparison\n');
fprintf('  relative RMS R             : %.6e\n', ...
    comparison.geometry.R_relative_rms);
fprintf('  absolute RMS Z             : %.6e\n', ...
    comparison.geometry.Z_absolute_rms);
fprintf('  relative RMS Bp            : %.6e\n', ...
    comparison.geometry.Bp_relative_rms);
fprintf('  Bp RMS on fitted boundary  : %.6e\n', ...
    comparison.geometry.Bp_fit_relative_rms);
fprintf('  relative RMS theta_psi_N   : %.6e\n', ...
    comparison.geometry.theta_psi_relative_rms);

fprintf('\nExtended coefficient comparison\n');
fprintf('  relative RMS g             : %.6e\n', ...
    comparison.coefficient.g_relative_rms);
fprintf('  relative RMS c             : %.6e\n', ...
    comparison.coefficient.c_relative_rms);
fprintf('  relative RMS f             : %.6e\n', ...
    comparison.coefficient.f_relative_rms);
fprintf('  relative RMS K_N           : %.6e\n', ...
    comparison.coefficient.K_relative_rms);

fprintf('\nCommon finite-element eigenproblem\n');
fprintf('  lambda direct metric       : %+.10e\n', lambda_metric);
fprintf('  lambda C-shape local       : %+.10e\n', lambda_local);
fprintf('  absolute / relative diff   : %.6e / %.6e\n', ...
    comparison.eigenvalue.absolute_difference, ...
    comparison.eigenvalue.relative_difference);
fprintf('  eigenfunction overlap      : %.8f\n', mode_overlap);

%% Comparison plots
figure('Color', 'w', 'Name', 'Direct metric vs C-shape geometry');
tiledlayout(2, 2, 'TileSpacing', 'compact');

nexttile;
plot(R_metric, Z_metric, '-k', 'LineWidth', 2);
hold on;
plot(R_local, Z_local, '--r', 'LineWidth', 1.5);
axis equal;
grid on;
xlabel('R/a_N');
ylabel('Z/a_N');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');
title('Target flux surface');

nexttile;
plot(theta_base, R_metric, '-k', 'LineWidth', 1.5);
hold on;
plot(theta_base, R_local, '--r', 'LineWidth', 1.5);
grid on;
xlabel('PEST angle \theta');
ylabel('R/a_N');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');

nexttile;
plot(theta_base, Bp_metric, '-k', 'LineWidth', 1.5);
hold on;
plot(theta_base, Bp_local, '--r', 'LineWidth', 1.5);
grid on;
xlabel('PEST angle \theta');
ylabel('B_p/B_N');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');

nexttile;
plot(theta_base, theta_psi_metric, '-k', 'LineWidth', 1.5);
hold on;
plot(theta_base, theta_psi_local_N, '--r', 'LineWidth', 1.5);
grid on;
xlabel('PEST angle \theta');
ylabel('\partial\theta/\partial\psi_N');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');

figure('Color', 'w', 'Name', 'Direct metric vs C-shape g-c-f');
tiledlayout(3, 1, 'TileSpacing', 'compact');

nexttile;
semilogy(theta_common, coeff_metric.g, '-k', 'LineWidth', 1.5);
hold on;
semilogy(theta_common, g_local, '--r', 'LineWidth', 1.5);
grid on;
ylabel('g');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');

nexttile;
plot(theta_common, coeff_metric.c, '-k', 'LineWidth', 1.5);
hold on;
plot(theta_common, c_local, '--r', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
ylabel('c');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');

nexttile;
semilogy(theta_common, coeff_metric.f, '-k', 'LineWidth', 1.5);
hold on;
semilogy(theta_common, f_local, '--r', 'LineWidth', 1.5);
grid on;
xlabel('Ballooning angle \theta');
ylabel('f');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');

figure('Color', 'w', 'Name', 'Direct metric vs C-shape eigenmode');
plot(theta_dof, mode_metric, '-k', 'LineWidth', 2);
hold on;
plot(theta_dof_local, mode_local, '--r', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
xlabel('Ballooning angle \theta');
ylabel('X/max|X|');
legend('Direct GS metric', 'C-shape local', 'Location', 'best');
title(sprintf(['Largest mode: \\lambda_{direct}=%.5g, ', ...
    '\\lambda_{Cshape}=%.5g'], lambda_metric, lambda_local));

%% Internal checks; differences between the two models are only reported
local_tolerance = 1.e-10;

assert(coeff_metric.diagnostics.B2MetricRelativeRms < 5.e-2, ...
    'The direct metric path has an inconsistent B^2 metric.');
assert(bal_local.validation.all_finite, ...
    'The C-shape local-equilibrium path produced NaN or Inf.');
assert(all(coeff_metric.g > 0) && all(coeff_metric.f > 0), ...
    'The direct metric path produced non-positive g or f.');
assert(all(g_local > 0) && all(f_local > 0), ...
    'The C-shape local-equilibrium path produced non-positive g or f.');
assert(all(isfinite([lambda_metric; lambda_local; mode_overlap])), ...
    'The common eigenvalue comparison produced a non-finite result.');
assert(abs(comparison.profile.q_metric - comparison.profile.q_local) ...
        < 1.e-8*max(1, abs(comparison.profile.q_metric)), ...
    'The two paths are not using the same q profile value.');
assert(comparison.local.minimum_flux_normal_projection > 0, ...
    'The C-shape Mercier normal is not aligned with grad(psi).');
assert(sign(merluc.dpsi_dr) == sign(param.dpsi_dr_fit), ...
    'The q-normalized and radial-fit dpsi/dr signs do not agree.');
assert(abs(comparison.local.theta_PEST_span - 2*pi) < local_tolerance, ...
    'The C-shape PEST angle does not span exactly 2*pi.');
assert(max([bal_local.validation.K_relative_error, ...
            bal_local.validation.D_relative_error, ...
            bal_local.validation.g_relative_error, ...
            bal_local.validation.c_relative_error, ...
            bal_local.validation.f_relative_error]) < local_tolerance, ...
    'A C-shape closed-form coefficient fails its vector identity.');

fprintf('\n  internal validation status : PASS\n');

function value = relative_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2)) ...
        /max(sqrt(mean(reference.^2)), eps);
end

function value = absolute_rms(actual, reference)
    actual = actual(:);
    reference = reference(:);
    value = sqrt(mean((actual - reference).^2));
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

function mode = normalize_mode(mode)
    [~, index] = max(abs(mode));
    mode = mode.*exp(-1i*angle(mode(index)));
    mode = real(mode);
    mode = mode/max(abs(mode));
end

function surf = reverse_closed_surface(surf)
    surf.R = [surf.R(1); flipud(surf.R(2:end))];
    surf.Z = [surf.Z(1); flipud(surf.Z(2:end))];

    R_next = circshift(surf.R, -1);
    Z_next = circshift(surf.Z, -1);
    surf.dl = hypot(R_next - surf.R, Z_next - surf.Z);
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
