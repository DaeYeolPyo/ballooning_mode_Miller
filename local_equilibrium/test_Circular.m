clc
close all
clearvars -except target_psiN

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));

%% Dedicated circular-equilibrium input
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_circular');
if ~exist('target_psiN', 'var') || isempty(target_psiN)
    target_psiN = 0.8;
end

eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, Npoints=1024);

%% Fit R = R0 + r*cos(theta), Z = r*sin(theta)
[param, bnd] = fit_Circular( ...
    eq, eqfunc, target_psiN, NTheta=1024);

figure('Color', 'w', 'Name', 'Circular surface fit');
plot(bnd.R, bnd.Z, '-k', 'LineWidth', 2);
hold on;
plot(bnd.Req, bnd.Zeq, '-r', 'LineWidth', 2);
plot(eq.rbbbs, eq.zbbbs, '-b', 'LineWidth', 1);
axis equal;
grid on;
xlabel('R [m]');
ylabel('Z [m]');
legend('Circular fit', 'GEQDSK target surface', 'Plasma boundary', ...
    'Location', 'best');
title(sprintf('Circular fit at \\psi_N = %.3f', target_psiN));

%% Exact shifted-circle Mercier-Luc field
circluc = Circular_Mercier_Luc(surf, param, bnd);
eq_on_circle = eqfunc.eval(circluc.R, circluc.Z);

Bp_geqdsk = eq_on_circle.Bp(:);
Bphi_geqdsk = eq_on_circle.Bphi(:);
B_geqdsk = sqrt(eq_on_circle.B2(:));
psiN_on_circle = eq_on_circle.psiN(:);

Bp_relerr = (circluc.Bp - Bp_geqdsk)./max(abs(Bp_geqdsk), eps);
Bphi_relerr = (circluc.Bphi - Bphi_geqdsk) ...
    ./max(abs(Bphi_geqdsk), eps);
B_relerr = (circluc.B - B_geqdsk)./max(abs(B_geqdsk), eps);

rms_Bp_relerr = sqrt(mean(Bp_relerr.^2));
max_Bp_relerr = max(abs(Bp_relerr));
rms_Bphi_relerr = sqrt(mean(Bphi_relerr.^2));
max_Bphi_relerr = max(abs(Bphi_relerr));
rms_B_relerr = sqrt(mean(B_relerr.^2));
max_B_relerr = max(abs(B_relerr));
psiN_fit_rms = sqrt(mean((psiN_on_circle - target_psiN).^2));
psiN_fit_max = max(abs(psiN_on_circle - target_psiN));

% In an exact shifted-circle equilibrium this inferred dpsi/dr is constant.
dpsi_dr_geqdsk = Bp_geqdsk.*circluc.R ...
    .*circluc.radial_metric_factor;
dpsi_dr_geqdsk_mean = mean(dpsi_dr_geqdsk);
dpsi_dr_geqdsk_nonuniformity = sqrt(mean(( ...
    dpsi_dr_geqdsk - dpsi_dr_geqdsk_mean).^2)) ...
    /max(abs(dpsi_dr_geqdsk_mean), eps);
dpsi_dr_mean_relerr = abs( ...
    circluc.dpsi_dr - dpsi_dr_geqdsk_mean) ...
    /max(abs(dpsi_dr_geqdsk_mean), eps);

figure('Color', 'w', 'Name', 'Circular vs GEQDSK poloidal field');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(circluc.theta_geo, circluc.Bp, '-k', 'LineWidth', 2);
hold on;
plot(circluc.theta_geo, Bp_geqdsk, '--r', 'LineWidth', 2);
grid on;
xlim([0, 2*pi]);
ylabel('B_p [T]');
legend('Exact shifted-circle', 'GEQDSK at circular (R,Z)', ...
    'Location', 'best');
title(sprintf('Poloidal field at \\psi_N = %.3f', target_psiN));
ylim([0.1 0.2])

nexttile;
plot(circluc.theta_geo, 100*Bp_relerr, '-b', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
xlim([0, 2*pi]);
xlabel('Circular angle \theta [rad]');
ylabel('Relative error [%]');

fprintf('Circular local-equilibrium fit at psiN = %.6f\n', target_psiN);
fprintf('  GEQDSK file                 : %s\n', geqdsk_file);
fprintf('  R0, r, A                    : %.8e, %.8e, %.8e\n', ...
    param.R0, param.r, param.A);
fprintf('  dR0/dr                      : %.8e\n', param.dR0_dr);
fprintf('  contour radial RMS / max   : %.6e / %.6e m\n', ...
    param.fit_residual_rms, param.fit_residual_max);
fprintf('  circular-path psiN range   : [%.8f, %.8f]\n', ...
    min(psiN_on_circle), max(psiN_on_circle));
fprintf('  RMS / max |psiN-psiN0|     : %.6e / %.6e\n', ...
    psiN_fit_rms, psiN_fit_max);
fprintf('  dpsi/dr exact circular     : %.8e Wb/m\n', circluc.dpsi_dr);
fprintf('  dpsi/dr GEQDSK mean        : %.8e Wb/m\n', ...
    dpsi_dr_geqdsk_mean);
fprintf('  dpsi/dr mean rel. error    : %.6e\n', ...
    dpsi_dr_mean_relerr);
fprintf('  dpsi/dr GEQDSK nonuniform. : %.6e\n', ...
    dpsi_dr_geqdsk_nonuniformity);
fprintf('  q target / q check         : %.10f / %.10f\n', ...
    surf.q, circluc.q_check);
fprintf('  analytic/numeric q-int err : %.6e\n', ...
    circluc.q_integral_relative_error);
fprintf('  closed-form Jacobian error : %.6e\n', ...
    circluc.jac_relative_error);
fprintf('  RMS / max Bp rel. error    : %.6e / %.6e\n', ...
    rms_Bp_relerr, max_Bp_relerr);
fprintf('  RMS / max Bphi rel. error  : %.6e / %.6e\n', ...
    rms_Bphi_relerr, max_Bphi_relerr);
fprintf('  RMS / max |B| rel. error   : %.6e / %.6e\n', ...
    rms_B_relerr, max_B_relerr);

%% Gamma and theta_psi on the circular local equilibrium
FFprime = eqfunc.FFprime(target_psiN);
pprime = eqfunc.pprime(target_psiN);
dq_dpsiN = eqfunc.qprimeN(target_psiN);
qprime = eqfunc.qprime(target_psiN);

[Gamma, theta_psi, theta_diag] = evaluate_theta_psi( ...
    circluc, FFprime, pprime, qprime, EnforcePeriodicity=true);

figure('Color', 'w', 'Name', 'Circular Gamma and theta-psi');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(circluc.theta_PEST, theta_diag.Gamma_input, '--r', ...
    'LineWidth', 1.5);
hold on;
plot(circluc.theta_PEST, Gamma, '-k', 'LineWidth', 2);
grid on;
xlim([0, 2*pi]);
ylabel('\Gamma [1/Wb]');
legend('Raw GEQDSK FF''', 'Circular-periodic FF''', ...
    'Location', 'best');

nexttile;
plot(circluc.theta_PEST, theta_diag.theta_psi_input, '--r', ...
    'LineWidth', 1.5);
hold on;
plot(circluc.theta_PEST, theta_psi, '-k', 'LineWidth', 2);
plot(circluc.theta_PEST, theta_diag.theta_psi_PEST, ':b', ...
    'LineWidth', 1.5);
grid on;
xlim([0, 2*pi]);
xlabel('PEST angle \theta [rad]');
ylabel('\theta_\psi [1/Wb]');
legend('Raw GEQDSK FF''', 'Periodic: geometric integral', ...
    'Periodic: PEST integral', 'Location', 'best');

theta_PEST_span = circluc.theta_PEST(end) - circluc.theta_PEST(1);
normal_projection_aligned = ...
    circluc.normal_sign*circluc.normal_projection_raw;

fprintf('\nCircular Gamma/theta_psi validation\n');
fprintf('  dq/dpsiN, dq/dpsi           : %.8e, %.8e 1/Wb\n', ...
    dq_dpsiN, qprime);
fprintf('  theta_PEST span / 2pi       : %.12f / %.12f\n', ...
    theta_PEST_span, 2*pi);
fprintf('  FFprime input / periodic    : %.8e / %.8e\n', ...
    theta_diag.FFprime_input, theta_diag.FFprime_periodic);
fprintf('  raw closure defect          : %.6e\n', ...
    theta_diag.input_closure_defect);
fprintf('  periodic closure (geo/PEST) : %+.8e / %+.8e 1/Wb\n', ...
    theta_diag.closure_used, theta_diag.closure_PEST);

%% Ballooning coefficients
theta_bal = 0.0;
theta_bnd = 5*pi;
a_N = param.r;
B_N = abs(param.B0);
bal_coef = evaluate_gcf(theta_diag, theta_bal, ...
    theta_bnd=theta_bnd, a_N=a_N, B_N=B_N);

figure('Color', 'w', 'Name', 'Circular ballooning coefficients');
tiledlayout(3, 1, 'TileSpacing', 'compact');
nexttile;
semilogy(bal_coef.theta, bal_coef.g, '-k', 'LineWidth', 1.5);
grid on;
ylabel('g');
nexttile;
plot(bal_coef.theta, bal_coef.c, '-r', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
ylabel('c');
nexttile;
semilogy(bal_coef.theta, bal_coef.f, '-b', 'LineWidth', 1.5);
grid on;
xlabel('Ballooning angle \theta');
ylabel('f');

gcf_validation = bal_coef.validation;
nperiod = numel(theta_diag.theta) - 1;
I_jump = bal_coef.I(1+nperiod:end) - bal_coef.I(1:end-nperiod);
I_jump_expected = 2*pi*qprime;
I_jump_error = max(abs(I_jump - I_jump_expected));
periodic_geometry_error = max([ ...
    max(abs(bal_coef.R(1+nperiod:end) ...
        - bal_coef.R(1:end-nperiod))), ...
    max(abs(bal_coef.Bp(1+nperiod:end) ...
        - bal_coef.Bp(1:end-nperiod))), ...
    max(abs(bal_coef.theta_psi(1+nperiod:end) ...
        - bal_coef.theta_psi(1:end-nperiod)))]);

zero_pressure_diag = theta_diag;
zero_pressure_diag.pprime = 0.0;
zero_pressure_bal = evaluate_gcf(zero_pressure_diag, theta_bal, ...
    theta_bnd=theta_bnd, a_N=a_N, B_N=B_N);
zero_pressure_c_error = max(abs(zero_pressure_bal.c));

fprintf('\nCircular ballooning g/c/f validation\n');
fprintf('  g min/max                   : %.8e, %.8e\n', ...
    min(bal_coef.g), max(bal_coef.g));
fprintf('  c min/max                   : %.8e, %.8e\n', ...
    min(bal_coef.c), max(bal_coef.c));
fprintf('  f min/max                   : %.8e, %.8e\n', ...
    min(bal_coef.f), max(bal_coef.f));
fprintf('  relative error K, D         : %.3e, %.3e\n', ...
    gcf_validation.K_relative_error, ...
    gcf_validation.D_relative_error);
fprintf('  relative error g, c, f      : %.3e, %.3e, %.3e\n', ...
    gcf_validation.g_relative_error, ...
    gcf_validation.c_relative_error, ...
    gcf_validation.f_relative_error);
fprintf('  I(theta+2pi) jump error     : %.3e\n', I_jump_error);
fprintf('  periodic geometry error     : %.3e\n', periodic_geometry_error);
fprintf('  pprime=0 max|c|             : %.3e\n', zero_pressure_c_error);

%% Machine-readable result and validation thresholds
circular_validation = struct();
circular_validation.geometry = struct( ...
    'radial_rms_m', param.fit_residual_rms, ...
    'radial_max_m', param.fit_residual_max, ...
    'psiN_rms', psiN_fit_rms, ...
    'psiN_max', psiN_fit_max);
circular_validation.field = struct( ...
    'Bp_relative_rms', rms_Bp_relerr, ...
    'Bp_relative_max', max_Bp_relerr, ...
    'Bphi_relative_rms', rms_Bphi_relerr, ...
    'B_relative_rms', rms_B_relerr, ...
    'dpsi_dr_nonuniformity', dpsi_dr_geqdsk_nonuniformity);
circular_validation.theta = struct( ...
    'closure_geo', theta_diag.closure_used, ...
    'closure_PEST', theta_diag.closure_PEST);
circular_validation.gcf = gcf_validation;

machine_tolerance = 1.e-10;
closure_tolerance = 1.e-10*max(1, max(abs(theta_psi)));
assert(circluc.jac_relative_error < machine_tolerance, ...
    'The circular closed-form Jacobian is inconsistent.');
assert(circluc.q_integral_relative_error < machine_tolerance, ...
    'The analytic and numerical circular q integrals disagree.');
assert(abs(circluc.q_check - surf.q) ...
        < machine_tolerance*max(1, abs(surf.q)), ...
    'The circular Bp does not reproduce the target q.');
assert(rms_Bp_relerr < 5.e-2, ...
    'Circular Bp validation failed: RMS relative error exceeds 5%%.');
assert(psiN_fit_rms < 2.e-2, ...
    'The GEQDSK target surface is not sufficiently circular.');
assert(all(normal_projection_aligned > 0), ...
    'The circular normal is not aligned with increasing r.');
assert(abs(theta_PEST_span - 2*pi) < machine_tolerance, ...
    'The circular PEST angle does not span exactly 2*pi.');
assert(abs(theta_diag.closure_used) < closure_tolerance, ...
    'The circular theta_psi is not periodic.');
assert(gcf_validation.all_finite ...
        && all(bal_coef.g > 0) && all(bal_coef.f > 0), ...
    'The circular ballooning coefficients are invalid.');
assert(max([gcf_validation.K_relative_error, ...
            gcf_validation.D_relative_error, ...
            gcf_validation.g_relative_error, ...
            gcf_validation.c_relative_error, ...
            gcf_validation.f_relative_error]) < machine_tolerance, ...
    'The circular g/c/f closed forms fail their vector checks.');
assert(I_jump_error < machine_tolerance*max(1, abs(I_jump_expected)), ...
    'The circular integrated shear has the wrong 2*pi jump.');
assert(periodic_geometry_error < machine_tolerance, ...
    'The circular periodic geometry does not repeat after 2*pi.');
assert(zero_pressure_c_error < machine_tolerance, ...
    'The circular pressure-curvature drive does not vanish for pprime=0.');

fprintf('  validation status           : PASS\n');
