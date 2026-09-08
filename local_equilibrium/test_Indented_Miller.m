clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(this_dir, '..', 'equilibrium'));

%% Dedicated indented-Miller equilibrium
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_I0.8');
target_psiN = 0.7;

eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);
surf = extract_flux_surface(eq, eqfunc, target_psiN, Npoints=1024);

%% Fit R0, a, delta, kappa, and I0
[param, bnd] = fit_Indented_Miller( ...
    eq, eqfunc, target_psiN, ...
    NTheta=256, NRadialFit=9, dpsi=3.e-2, ...
    MaxIterations=250, NPhaseStarts=32, NXStarts=33);
%%

figure('Color', 'w', 'Name', 'Indented-Miller surface fit');
plot(bnd.R, bnd.Z, '-k', 'LineWidth', 2);
hold on;
plot(bnd.Req, bnd.Zeq, '--r', 'LineWidth', 1.5);
plot(eq.rbbbs, eq.zbbbs, '-b', 'LineWidth', 1);
grid on;
xlabel('R [m]');
ylabel('Z [m]');
legend('Indented-Miller fit', 'GEQDSK target surface', ...
    'Plasma boundary', 'Location', 'best');

%% Independent flux-surface check on the fitted analytic curve
eq_on_fit = eqfunc.eval(bnd.R, bnd.Z);
psiN_error = eq_on_fit.psiN(:) - target_psiN;
psiN_fit_rms = sqrt(mean(psiN_error.^2));
psiN_fit_max = max(abs(psiN_error));

figure('Color', 'w', 'Name', 'Indented-Miller fit residuals');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(bnd.theta, eq_on_fit.psiN, '-k', 'LineWidth', 1.5);
hold on;
yline(target_psiN, '--r');
grid on;
xlim([0, 2*pi]);
ylabel('\psi_N(R(\theta),Z(\theta))');
title('Flux value sampled on the analytic fitted curve');

nexttile;
plot(bnd.theta, 1.e3*psiN_error, '-b', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
xlim([0, 2*pi]);
xlabel('Analytic angle \theta [rad]');
ylabel('10^3(\psi_N-\psi_{N,target})');

%% Check that neighboring surfaces give smooth radial shape parameters
radial = param.radial_fit;
figure('Color', 'w', 'Name', 'Indented-Miller radial fit');
tiledlayout(2, 2, 'TileSpacing', 'compact');

nexttile;
plot(radial.psiN, radial.R0, 'ok-', 'LineWidth', 1.2, ...
    'MarkerFaceColor', 'w');
grid on;
xlabel('\psi_N');
ylabel('R_0 [m]');

nexttile;
plot(radial.psiN, radial.delta, 'ok-', 'LineWidth', 1.2, ...
    'MarkerFaceColor', 'w');
hold on;
yline(-0.6, ':', 'Nominal LCFS value');
grid on;
xlabel('\psi_N');
ylabel('\delta');

nexttile;
plot(radial.psiN, radial.kappa, 'ok-', 'LineWidth', 1.2, ...
    'MarkerFaceColor', 'w');
hold on;
yline(1.5, ':', 'Nominal LCFS value');
grid on;
xlabel('\psi_N');
ylabel('\kappa');

nexttile;
plot(radial.psiN, radial.I0, 'ok-', 'LineWidth', 1.2, ...
    'MarkerFaceColor', 'w');
hold on;
yline(0.8, ':', 'Nominal LCFS value');
grid on;
xlabel('\psi_N');
ylabel('I_0');

%% Report and validate
radial_converged = cellfun(@(item) item.converged, radial.diagnostics);
radial_relative_rms = cellfun( ...
    @(item) item.relative_rms_distance, radial.diagnostics);

fprintf('Indented-Miller fit: %s\n', geqdsk_file);
fprintf('  target psiN                 : %.6f\n', target_psiN);
fprintf('  R0, a                       : %.10f, %.10f m\n', ...
    param.R0, param.a);
fprintf('  x=asin(delta), b, d         : %.10f, %.10f, %.10f\n', ...
    param.x, param.b, param.d);
fprintf('  delta, kappa, I0            : %.10f, %.10f, %.10f\n', ...
    param.delta, param.kappa, param.I0);
fprintf('  geometric RMS / max [m]     : %.6e / %.6e\n', ...
    param.fit.rms_distance, param.fit.max_distance);
fprintf('  relative geometric RMS      : %.6e\n', ...
    param.fit.relative_rms_distance);
fprintf('  fitted-curve psiN RMS / max : %.6e / %.6e\n', ...
    psiN_fit_rms, psiN_fit_max);
fprintf('  center fit iterations       : %d (converged=%d)\n', ...
    param.fit.iterations, param.fit.converged);
fprintf('  radial fits converged       : %d / %d\n', ...
    nnz(radial_converged), numel(radial_converged));
fprintf('  max radial relative RMS     : %.6e\n', ...
    max(radial_relative_rms));
fprintf('  dR0/da, dx/da               : %.8e, %.8e 1/m\n', ...
    param.dR0_da, param.dx_da);
fprintf('  db/da, dd/da                : %.8e, %.8e\n', ...
    param.db_da, param.dd_da);
fprintf('  s_delta, s_kappa, s_I0      : %.8e, %.8e, %.8e\n', ...
    param.s_delta, param.s_kappa, param.s_I0);
fprintf('  dpsi/da from radial fit      : %.8e Wb/m\n', ...
    param.dpsi_da_fit);

assert(param.a > 0 && param.kappa > 0 && param.I0 >= 0 ...
    && abs(param.delta) < 1, ...
    'The fitted indented-Miller parameters are not physical.');
assert(all(radial_converged), ...
    'At least one neighboring-surface fit did not converge.');
assert(all(diff(radial.a) > 0), ...
    'The fitted minor radius is not monotone with psiN.');
assert(param.fit.relative_rms_distance < 2.e-3, ...
    'The target-surface relative geometric RMS error exceeds 0.2%%.');
assert(max(radial_relative_rms) < 3.e-3, ...
    'A radial surface relative geometric RMS error exceeds 0.3%%.');
assert(psiN_fit_rms < 2.e-3 && psiN_fit_max < 5.e-3, ...
    'The fitted curve does not remain sufficiently close to target psiN.');

indent_fit_validation = struct();
indent_fit_validation.geqdsk_file = geqdsk_file;
indent_fit_validation.target_psiN = target_psiN;
indent_fit_validation.param = param;
indent_fit_validation.bnd = bnd;
indent_fit_validation.surf = surf;
indent_fit_validation.psiN_fit_rms = psiN_fit_rms;
indent_fit_validation.psiN_fit_max = psiN_fit_max;
indent_fit_validation.radial_converged = radial_converged;
indent_fit_validation.radial_relative_rms = radial_relative_rms;
indent_fit_validation.fitting_status = 'PASS';

fprintf('  fitting validation status   : PASS\n');

%% Evaluate the Mercier-Luc poloidal field
merluc = Indented_Miller_Mercier_Luc(surf, param, bnd);
eq_on_shape = eqfunc.eval(merluc.R, merluc.Z);

Bp_geqdsk = eq_on_shape.Bp(:);
Bphi_geqdsk = eq_on_shape.Bphi(:);
B_geqdsk = sqrt(eq_on_shape.B2(:));
psiN_on_shape = eq_on_shape.psiN(:);

Bp_relerr = (merluc.Bp - Bp_geqdsk)./max(abs(Bp_geqdsk), eps);
Bphi_relerr = (merluc.Bphi - Bphi_geqdsk) ...
    ./max(abs(Bphi_geqdsk), eps);
B_relerr = (merluc.B - B_geqdsk)./max(abs(B_geqdsk), eps);

rms_Bp_relerr = sqrt(mean(Bp_relerr.^2));
max_Bp_relerr = max(abs(Bp_relerr));
rms_Bphi_relerr = sqrt(mean(Bphi_relerr.^2));
max_Bphi_relerr = max(abs(Bphi_relerr));
rms_B_relerr = sqrt(mean(B_relerr.^2));
max_B_relerr = max(abs(B_relerr));
dpsi_da_relerr = abs(merluc.dpsi_da - param.dpsi_da_fit) ...
    /max(abs(param.dpsi_da_fit), eps);

figure('Color', 'w', 'Name', ...
    'Indented-Miller vs GEQDSK magnetic field');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(merluc.theta_geo, merluc.Bp, '-k', 'LineWidth', 2);
hold on;
plot(merluc.theta_geo, Bp_geqdsk, '--r', 'LineWidth', 2);
grid on;
xlim([0, 2*pi]);
ylabel('B_p [T]');
legend('Indented-Miller/Mercier-Luc', ...
    'GEQDSK at fitted (R,Z)', 'Location', 'best');
title(sprintf('Poloidal field at \\psi_N = %.3f', target_psiN));

nexttile;
plot(merluc.theta_geo, 100*Bp_relerr, '-b', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
xlim([0, 2*pi]);
xlabel('Indented-Miller angle \theta [rad]');
ylabel('Relative error [%]');

fprintf('\nIndented-Miller vs GEQDSK magnetic field\n');
fprintf('  dpsi/da from q constraint   : %.8e Wb/m\n', ...
    merluc.dpsi_da);
fprintf('  dpsi/da from radial fit     : %.8e Wb/m\n', ...
    param.dpsi_da_fit);
fprintf('  dpsi/da relative difference : %.6e\n', dpsi_da_relerr);
fprintf('  q target / q check          : %.10f / %.10f\n', ...
    surf.q, merluc.q_check);
fprintf('  fitted-surface psiN range   : [%.8f, %.8f]\n', ...
    min(psiN_on_shape), max(psiN_on_shape));
fprintf('  RMS / max Bp rel. error     : %.6e / %.6e\n', ...
    rms_Bp_relerr, max_Bp_relerr);
fprintf('  RMS / max Bphi rel. error   : %.6e / %.6e\n', ...
    rms_Bphi_relerr, max_Bphi_relerr);
fprintf('  RMS / max |B| rel. error    : %.6e / %.6e\n', ...
    rms_B_relerr, max_B_relerr);

%% Evaluate Gamma and theta_psi
FFprime = eqfunc.FFprime(target_psiN);
pprime = eqfunc.pprime(target_psiN);
dpsiN = 0.01;
dq_dpsiN = (eqfunc.q(target_psiN + dpsiN) ...
    - eqfunc.q(target_psiN - dpsiN))/(2*dpsiN);
qprime = dq_dpsiN/(eq.sibry - eq.simag);

[Gamma, theta_psi, theta_diag] = evaluate_theta_psi( ...
    merluc, FFprime, pprime, qprime, EnforcePeriodicity=true);

figure('Color', 'w', 'Name', ...
    'Indented-Miller Gamma and theta-psi periodicity');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(merluc.theta_PEST, theta_diag.Gamma_input, '--r', ...
    'LineWidth', 1.5);
hold on;
plot(merluc.theta_PEST, Gamma, '-k', 'LineWidth', 2);
grid on;
xlim([0, 2*pi]);
ylabel('\Gamma [1/Wb]');
legend('Raw GEQDSK FF''', 'Indented-Miller-periodic FF''', ...
    'Location', 'best');
title('Radial derivative of the PEST-angle mapping');

nexttile;
plot(merluc.theta_PEST, theta_diag.theta_psi_input, '--r', ...
    'LineWidth', 1.5);
hold on;
plot(merluc.theta_PEST, theta_psi, '-k', 'LineWidth', 2);
plot(merluc.theta_PEST, theta_diag.theta_psi_PEST, ':b', ...
    'LineWidth', 1.5);
grid on;
xlim([0, 2*pi]);
xlabel('PEST angle \theta [rad]');
ylabel('\theta_\psi [1/Wb]');
legend('Raw GEQDSK FF''', 'Periodic: geometric-angle integral', ...
    'Periodic: PEST-angle integral', 'Location', 'best');

theta_PEST_span = merluc.theta_PEST(end) - merluc.theta_PEST(1);
normal_projection_r = ...
    merluc.normal_sign*merluc.normal_projection_raw;
normal_projection_psi = sign(merluc.dpsi_da)*normal_projection_r;

fprintf('\nGamma/theta_psi periodicity validation\n');
fprintf('  psi scale                   : %.8e Wb/rad\n', ...
    eq.sibry - eq.simag);
fprintf('  dq/dpsiN, dq/dpsi           : %.8e, %.8e 1/Wb\n', ...
    dq_dpsiN, qprime);
fprintf('  Mercier normal sign         : %+d\n', merluc.normal_sign);
fprintf('  theta_PEST span / 2pi       : %.12f / %.12f\n', ...
    theta_PEST_span, 2*pi);
fprintf('  FFprime GEQDSK              : %.8e\n', ...
    theta_diag.FFprime_input);
fprintf('  FFprime periodic            : %.8e\n', ...
    theta_diag.FFprime_periodic);
fprintf('  raw closure                 : %+.8e 1/Wb\n', ...
    theta_diag.closure_input);
fprintf('  raw closure defect          : %.6e\n', ...
    theta_diag.input_closure_defect);
fprintf('  periodic closure (geo/PEST) : %+.3e / %+.3e 1/Wb\n', ...
    theta_diag.closure_used, theta_diag.closure_PEST);
fprintf('  geo/PEST max difference     : %.8e 1/Wb\n', ...
    theta_diag.integration_max_abs_difference);

closure_tolerance = 1.e-10*max(1, max(abs(theta_psi)));
assert(rms_Bp_relerr < 2.e-2, ...
    'Indented-Miller Bp RMS relative error exceeds 2%%.');
assert(rms_Bphi_relerr < 1.e-2 && rms_B_relerr < 1.e-2, ...
    'Indented-Miller toroidal/total-field RMS error exceeds 1%%.');
assert(dpsi_da_relerr < 2.e-2, ...
    'The q-constrained and radial-fit dpsi/da differ by more than 2%%.');
assert(all(normal_projection_psi > 0), ...
    'The Mercier normal is not aligned with grad(psi).');
assert(abs(theta_PEST_span - 2*pi) < 1.e-10, ...
    'The PEST angle does not span exactly 2*pi.');
assert(abs(theta_diag.closure_used) < closure_tolerance, ...
    'Indented-Miller-consistent theta_psi is not periodic.');

fprintf('  field/theta validation      : PASS\n');

%% Evaluate ballooning-equation coefficients g, c, and f
theta_bal = 0.0;
theta_bnd = 5*pi;
a_N = param.a;
B_N = abs(param.B0);

bal_coef = evaluate_gcf(theta_diag, theta_bal, ...
    theta_bnd=theta_bnd, a_N=a_N, B_N=B_N);

figure('Color', 'w', 'Name', ...
    'Indented-Miller ballooning coefficients');
tiledlayout(2, 2, 'TileSpacing', 'compact');

nexttile;
semilogy(bal_coef.theta, bal_coef.g, '-k', 'LineWidth', 1.5);
hold on;
semilogy(bal_coef.theta, bal_coef.f, '--b', 'LineWidth', 1.5);
grid on;
xlabel('Ballooning angle \theta');
ylabel('Positive coefficient');
legend('g', 'f', 'Location', 'best');
title('Field-line bending and inertia');

nexttile;
plot(bal_coef.theta, bal_coef.c, '-r', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
xlabel('Ballooning angle \theta');
ylabel('c');
title('Pressure-curvature drive');

nexttile;
plot(bal_coef.theta, bal_coef.I, '-m', 'LineWidth', 1.5);
grid on;
xlabel('Ballooning angle \theta');
ylabel('I');
title('Integrated local shear');

nexttile;
yyaxis left;
semilogy(bal_coef.theta, bal_coef.K, '-k', 'LineWidth', 1.5);
ylabel('K = |\nabla\alpha|^2');
yyaxis right;
plot(bal_coef.theta, bal_coef.D, '-g', 'LineWidth', 1.2);
ylabel('D');
grid on;
xlabel('Ballooning angle \theta');
title('Metric and curvature contraction');

gcf_validation = bal_coef.validation;
nperiod = numel(theta_diag.theta) - 1;
I_jump = bal_coef.I(1+nperiod:end) - bal_coef.I(1:end-nperiod);
I_jump_expected = 2*pi*qprime;
I_jump_error = max(abs(I_jump - I_jump_expected));

R_periodic_error = max(abs( ...
    bal_coef.R(1+nperiod:end) - bal_coef.R(1:end-nperiod)));
Bp_periodic_error = max(abs( ...
    bal_coef.Bp(1+nperiod:end) - bal_coef.Bp(1:end-nperiod)));
theta_psi_periodic_error = max(abs( ...
    bal_coef.theta_psi(1+nperiod:end) ...
    - bal_coef.theta_psi(1:end-nperiod)));

zero_pressure_diag = theta_diag;
zero_pressure_diag.pprime = 0.0;
zero_pressure_bal = evaluate_gcf(zero_pressure_diag, theta_bal, ...
    theta_bnd=theta_bnd, a_N=a_N, B_N=B_N);
zero_pressure_c_error = max(abs(zero_pressure_bal.c));

fprintf('\nBallooning g/c/f validation\n');
fprintf('  theta range                 : [%.8f, %.8f]\n', ...
    bal_coef.theta(1), bal_coef.theta(end));
fprintf('  a_N, B_N                    : %.8e m, %.8e T\n', ...
    a_N, B_N);
fprintf('  d(q)/d(psi_N,Gaur)          : %.8e\n', ...
    bal_coef.dq_dpsiN);
fprintf('  d(mu0*p/B_N^2)/dpsi_N       : %.8e\n', ...
    bal_coef.dpbar_dpsiN);
fprintf('  g min/max                   : %.8e, %.8e\n', ...
    min(bal_coef.g), max(bal_coef.g));
fprintf('  c min/max                   : %.8e, %.8e\n', ...
    min(bal_coef.c), max(bal_coef.c));
fprintf('  f min/max                   : %.8e, %.8e\n', ...
    min(bal_coef.f), max(bal_coef.f));
fprintf('  basis norm/orth errors      : %.3e, %.3e, %.3e\n', ...
    gcf_validation.basis_normal_norm_error, ...
    gcf_validation.basis_tangent_norm_error, ...
    gcf_validation.basis_orthogonality_error);
fprintf('  relative error K, D         : %.3e, %.3e\n', ...
    gcf_validation.K_relative_error, ...
    gcf_validation.D_relative_error);
fprintf('  relative error g, c, f      : %.3e, %.3e, %.3e\n', ...
    gcf_validation.g_relative_error, ...
    gcf_validation.c_relative_error, ...
    gcf_validation.f_relative_error);
fprintf('  I(theta+2pi) jump error     : %.3e\n', I_jump_error);
fprintf('  periodic R/Bp/thetapsi err  : %.3e, %.3e, %.3e\n', ...
    R_periodic_error, Bp_periodic_error, theta_psi_periodic_error);
fprintf('  pprime=0 max|c|             : %.3e\n', ...
    zero_pressure_c_error);

gcf_tolerance = 1.e-10;
periodic_tolerance = 1.e-10;
assert(gcf_validation.all_finite, ...
    'The ballooning coefficients contain NaN or Inf.');
assert(all(bal_coef.g > 0) && all(bal_coef.f > 0), ...
    'Expected positive g and f for F/q > 0.');
assert(gcf_validation.K_relative_error < gcf_tolerance, ...
    'Closed-form K does not match its vector definition.');
assert(gcf_validation.D_relative_error < gcf_tolerance, ...
    'Closed-form D does not match its vector definition.');
assert(gcf_validation.g_relative_error < gcf_tolerance ...
        && gcf_validation.c_relative_error < gcf_tolerance ...
        && gcf_validation.f_relative_error < gcf_tolerance, ...
    'Closed-form g/c/f do not match the vector definitions.');
assert(I_jump_error < periodic_tolerance*max(1, abs(I_jump_expected)), ...
    'Integrated local shear does not have the expected 2*pi*qprime jump.');
assert(max([R_periodic_error, Bp_periodic_error, ...
        theta_psi_periodic_error]) < periodic_tolerance, ...
    'A periodic geometry field does not repeat after 2*pi.');
assert(zero_pressure_c_error < gcf_tolerance, ...
    'The coefficient c does not vanish when pprime=0.');

indent_fit_validation.merluc = merluc;
indent_fit_validation.theta_diag = theta_diag;
indent_fit_validation.bal_coef = bal_coef;
indent_fit_validation.rms_Bp_relerr = rms_Bp_relerr;
indent_fit_validation.max_Bp_relerr = max_Bp_relerr;
indent_fit_validation.dpsi_da_relerr = dpsi_da_relerr;
indent_fit_validation.gcf_validation = gcf_validation;
indent_fit_validation.status = 'PASS';

fprintf('  ballooning validation       : PASS\n');
fprintf('\nOverall validation status     : PASS\n');
