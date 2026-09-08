function miller_validation = test_Miller( ...
        geqdsk_file, target_psiN, miller_fit_dpsi)
    clc
    close all

    this_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(this_dir, '..', 'equilibrium'));

%% Load equilibrium from GEQDSK
if nargin < 1 || isempty(geqdsk_file)
    geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_PT0.6');
end
eq = read_geqdsk(geqdsk_file);
eqfunc = build_interpolants(eq);

%% Set target flux surface
if nargin < 2 || isempty(target_psiN)
    target_psiN = 0.7;
end
if nargin < 3 || isempty(miller_fit_dpsi)
    miller_fit_dpsi = 3.e-2;
end

surf = extract_flux_surface(eq, eqfunc, target_psiN, Npoints=1024);

%% Fit the flux surface into analytic shape
[param, bnd] = fit_Miller(eq, eqfunc, target_psiN, ...
    NTheta=1024, dpsi=miller_fit_dpsi);

figure('Color', 'w', 'Name', 'Miller surface fit');
plot(bnd.R, bnd.Z, '-k', 'LineWidth', 2);
hold on;
plot(bnd.Req, bnd.Zeq, '-r', 'LineWidth', 2);
plot(eq.rbbbs, eq.zbbbs, '-b', 'LineWidth', 1);
grid on;
xlabel('R [m]');
ylabel('Z [m]');
legend('Miller fit', 'GEQDSK target surface', 'Plasma boundary', ...
    'Location', 'best');

%% Evaluate poloidal field
merluc = Miller_Mercier_Luc(surf, param, bnd);

% Evaluate the GEQDSK field at the same (R,Z) points as the Miller fit.
eq_on_Miller = eqfunc.eval(merluc.R, merluc.Z);
Bp_geqdsk = eq_on_Miller.Bp(:);
Bphi_geqdsk = eq_on_Miller.Bphi(:);
B_geqdsk = sqrt(eq_on_Miller.B2(:));
psiN_on_Miller = eq_on_Miller.psiN(:);

relerr = (merluc.Bp - Bp_geqdsk)./max(abs(Bp_geqdsk), eps);
rms_Bp_relerr = sqrt(mean(relerr.^2));
max_Bp_relerr = max(abs(relerr));
relerr_Bphi = (merluc.Bphi - Bphi_geqdsk) ...
    ./max(abs(Bphi_geqdsk), eps);
rms_Bphi_relerr = sqrt(mean(relerr_Bphi.^2));
max_Bphi_relerr = max(abs(relerr_Bphi));
relerr_B = (merluc.B - B_geqdsk)./max(abs(B_geqdsk), eps);
rms_B_relerr = sqrt(mean(relerr_B.^2));
max_B_relerr = max(abs(relerr_B));
psiN_fit_rms = sqrt(mean((psiN_on_Miller - target_psiN).^2));
psiN_fit_max = max(abs(psiN_on_Miller - target_psiN));

figure('Color', 'w', 'Name', 'Miller vs GEQDSK poloidal field');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(merluc.theta_geo, merluc.Bp, '-k', 'LineWidth', 2);
hold on;
plot(merluc.theta_geo, Bp_geqdsk, '--r', 'LineWidth', 2);
grid on;
xlim([0, 2*pi]);
ylabel('B_p [T]');
legend('Miller/Mercier-Luc', 'GEQDSK at Miller (R,Z)', ...
    'Location', 'best');
title(sprintf('Poloidal field at \\psi_N = %.3f', target_psiN));

nexttile;
plot(merluc.theta_geo, 100*relerr, '-b', 'LineWidth', 1.5);
yline(0, ':k');
grid on;
xlim([0, 2*pi]);
xlabel('Miller geometric angle \theta [rad]');
ylabel('Relative error [%]');

fprintf('Miller vs GEQDSK Bp at psiN = %.6f\n', target_psiN);
fprintf('  GEQDSK file               : %s\n', geqdsk_file);
fprintf('  R0, r, kappa, delta       : %.8e, %.8e, %.8e, %.8e\n', ...
    param.R0, param.r, param.kappa, param.delta);
fprintf('  dpsi/dr from q constraint : %.8e Wb/m\n', merluc.dpsi_dr);
fprintf('  q target / q check        : %.8f / %.8f\n', ...
    surf.q, merluc.q_check);
fprintf('  fitted-surface psiN range : [%.8f, %.8f]\n', ...
    min(psiN_on_Miller), max(psiN_on_Miller));
fprintf('  RMS / max |psiN-psiN0|    : %.6e / %.6e\n', ...
    psiN_fit_rms, psiN_fit_max);
fprintf('  RMS relative Bp error     : %.6e\n', ...
    rms_Bp_relerr);
fprintf('  Max relative Bp error     : %.6e\n', ...
    max_Bp_relerr);
fprintf('  RMS / max Bphi rel. error : %.6e / %.6e\n', ...
    rms_Bphi_relerr, max_Bphi_relerr);
fprintf('  RMS / max |B| rel. error  : %.6e / %.6e\n', ...
    rms_B_relerr, max_B_relerr);
fprintf('  radial-fit half width     : %.6f in psiN\n', ...
    param.radial_fit.half_width);
fprintf('  radial-fit surfaces       : %d\n', ...
    param.radial_fit.npoints);
fprintf('  dR0/dr, s_kappa, s_delta  : %.8f, %.8f, %.8f\n', ...
    param.dR0_dr, param.s_kappa, param.s_delta);

%% Evaluate Gamma and theta_psi
FFprime = eqfunc.FFprime(target_psiN);
pprime = eqfunc.pprime(target_psiN);

dpsiN = 0.01;
dq_dpsiN = (eqfunc.q(target_psiN + dpsiN) ...
    - eqfunc.q(target_psiN - dpsiN))/(2*dpsiN);
qprime = dq_dpsiN/(eq.sibry - eq.simag);
[Gamma, theta_psi, theta_diag] = evaluate_theta_psi( ...
    merluc, FFprime, pprime, qprime, EnforcePeriodicity=true);

figure('Color', 'w', 'Name', 'Gamma and theta-psi periodicity');
tiledlayout(2, 1, 'TileSpacing', 'compact');

nexttile;
plot(merluc.theta_PEST, theta_diag.Gamma_input, '--r', ...
    'LineWidth', 1.5);
hold on;
plot(merluc.theta_PEST, Gamma, '-k', 'LineWidth', 2);
grid on;
xlim([0, 2*pi]);
ylabel('\Gamma [1/Wb]');
legend('Raw GEQDSK FF''', 'Miller-periodic FF''', ...
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
normal_projection_aligned = ...
    merluc.normal_sign*merluc.normal_projection_raw;

fprintf('\nGamma/theta_psi periodicity validation\n');
fprintf('  psi scale                 : %.8e Wb/rad\n', ...
    eq.sibry - eq.simag);
fprintf('  dq/dpsiN, dq/dpsi         : %.8e, %.8e 1/Wb\n', ...
    dq_dpsiN, qprime);
fprintf('  Mercier normal sign       : %+d\n', merluc.normal_sign);
fprintf('  theta_PEST span / 2pi     : %.12f / %.12f\n', ...
    theta_PEST_span, 2*pi);
fprintf('  FFprime GEQDSK            : %.8e\n', ...
    theta_diag.FFprime_input);
fprintf('  FFprime Miller-periodic   : %.8e\n', ...
    theta_diag.FFprime_periodic);
fprintf('  raw closure               : %+.8e 1/Wb\n', ...
    theta_diag.closure_input);
fprintf('  raw closure defect        : %.6e\n', ...
    theta_diag.input_closure_defect);
fprintf('  periodic closure (geo)    : %+.8e 1/Wb\n', ...
    theta_diag.closure_used);
fprintf('  periodic closure (PEST)   : %+.8e 1/Wb\n', ...
    theta_diag.closure_PEST);
fprintf('  geo/PEST max difference   : %.8e 1/Wb\n', ...
    theta_diag.integration_max_abs_difference);

closure_tolerance = 1.e-10*max(1, max(abs(theta_psi)));
assert(rms_Bp_relerr < 5.e-2, ...
    'Miller Bp validation failed: RMS relative error exceeds 5%%.');
assert(rms_Bphi_relerr < 5.e-2 && rms_B_relerr < 5.e-2, ...
    'Miller toroidal/total-field validation failed: RMS error exceeds 5%%.');
assert(all(normal_projection_aligned > 0), ...
    'Mercier normal is not aligned with increasing Miller r.');
assert(abs(theta_PEST_span - 2*pi) < 1.e-10, ...
    'PEST angle does not span exactly 2*pi.');
assert(abs(theta_diag.closure_used) < closure_tolerance, ...
    'Miller-consistent theta_psi is not periodic.');

fprintf('  validation status         : PASS\n');

%% Evaluate ballooning coefficients g, c, f
theta_bal = 0.0;
theta_bnd = 5*pi;
a_N = param.r;
B_N = abs(param.B0);

bal_coef = evaluate_gcf(theta_diag, theta_bal, ...
    theta_bnd=theta_bnd, a_N=a_N, B_N=B_N);

figure('Color', 'w', 'Name', 'Miller ballooning coefficients');
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
fprintf('  theta range               : [%.8f, %.8f]\n', ...
    bal_coef.theta(1), bal_coef.theta(end));
fprintf('  a_N, B_N                  : %.8e m, %.8e T\n', a_N, B_N);
fprintf('  d(q)/d(psi_N,Gaur)        : %.8e\n', bal_coef.dq_dpsiN);
fprintf('  d(mu0*p/B_N^2)/dpsi_N    : %.8e\n', bal_coef.dpbar_dpsiN);
fprintf('  g min/max                 : %.8e, %.8e\n', ...
    min(bal_coef.g), max(bal_coef.g));
fprintf('  c min/max                 : %.8e, %.8e\n', ...
    min(bal_coef.c), max(bal_coef.c));
fprintf('  f min/max                 : %.8e, %.8e\n', ...
    min(bal_coef.f), max(bal_coef.f));
fprintf('  basis norm/orth errors    : %.3e, %.3e, %.3e\n', ...
    gcf_validation.basis_normal_norm_error, ...
    gcf_validation.basis_tangent_norm_error, ...
    gcf_validation.basis_orthogonality_error);
fprintf('  relative error K, D       : %.3e, %.3e\n', ...
    gcf_validation.K_relative_error, ...
    gcf_validation.D_relative_error);
fprintf('  relative error g, c, f    : %.3e, %.3e, %.3e\n', ...
    gcf_validation.g_relative_error, ...
    gcf_validation.c_relative_error, ...
    gcf_validation.f_relative_error);
fprintf('  I(theta+2pi) jump error   : %.3e\n', I_jump_error);
fprintf('  periodic R/Bp/thetapsi err: %.3e, %.3e, %.3e\n', ...
    R_periodic_error, Bp_periodic_error, theta_psi_periodic_error);
fprintf('  pprime=0 max|c|           : %.3e\n', zero_pressure_c_error);

gcf_tolerance = 1.e-10;
periodic_tolerance = 1.e-10;
assert(gcf_validation.all_finite, ...
    'The ballooning coefficients contain NaN or Inf.');
assert(all(bal_coef.g > 0) && all(bal_coef.f > 0), ...
    'Expected positive g and f for F/q > 0.');
assert(gcf_validation.K_relative_error < gcf_tolerance, ...
    'Closed-form K does not match the vector definition.');
assert(gcf_validation.D_relative_error < gcf_tolerance, ...
    'Closed-form D does not match the vector cross product.');
assert(gcf_validation.g_relative_error < gcf_tolerance ...
        && gcf_validation.c_relative_error < gcf_tolerance ...
        && gcf_validation.f_relative_error < gcf_tolerance, ...
    'Closed-form g/c/f do not match the Gaur definitions.');
assert(I_jump_error < periodic_tolerance*max(1, abs(I_jump_expected)), ...
    'Integrated local shear does not have the expected 2*pi*qprime jump.');
assert(max([R_periodic_error, Bp_periodic_error, ...
        theta_psi_periodic_error]) < periodic_tolerance, ...
    'A periodic geometry field does not repeat after 2*pi.');
assert(zero_pressure_c_error < gcf_tolerance, ...
    'The pressure-curvature coefficient c does not vanish for pprime=0.');

fprintf('  validation status         : PASS\n');

miller_validation = struct();
miller_validation.model = 'miller-mercier-luc';
miller_validation.geqdsk_file = geqdsk_file;
miller_validation.target_psiN = target_psiN;
miller_validation.miller_fit_dpsi = miller_fit_dpsi;
miller_validation.param = param;
miller_validation.bnd = bnd;
miller_validation.surf = surf;
miller_validation.merluc = merluc;
miller_validation.theta_diag = theta_diag;
miller_validation.bal_coef = bal_coef;
miller_validation.rms_Bp_relerr = rms_Bp_relerr;
miller_validation.max_Bp_relerr = max_Bp_relerr;
miller_validation.psiN_fit_rms = psiN_fit_rms;
miller_validation.psiN_fit_max = psiN_fit_max;
miller_validation.status = 'PASS';
end
