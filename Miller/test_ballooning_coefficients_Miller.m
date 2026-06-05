clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'salpha'));

p = DshapeMillerParams();

theta0 = 0.0;
nTheta = 257;
nGeom = 1201;

bal = miller_ballooning_coefficients(p, ...
    'Theta0', theta0, ...
    'NTheta', nTheta, ...
    'NGeom', nGeom);

theta = bal.theta(:);

% Recompute g, c, f from the diagnostic fields stored in bal.  This checks
% that the returned coefficients are consistent with the returned geometry.
thetaShift = theta - bal.theta0;
gradNAlphaCheck = bal.gradN_zeta ...
    - bal.q.*bal.gradN_theta ...
    - bal.dq_dpsin.*thetaShift.*bal.gradN_psin;
gradNAlphaSqCheck = sum(gradNAlphaCheck.^2, 2);
crossTermCheck = sum(cross(bal.bvec, bal.gradN_H, 2).*gradNAlphaCheck, 2);

gCheck = bal.b_dot_gradN_theta.*gradNAlphaSqCheck./bal.BoverBN;
fCheck = gradNAlphaSqCheck./(bal.b_dot_gradN_theta.*bal.BoverBN.^3);
cCheck = 2./bal.b_dot_gradN_theta./bal.BoverBN.^4 ...
    .*bal.dpbar_dpsin.*crossTermCheck;

rel = @(a,b) max(abs(a(:)-b(:)), [], 'omitnan') ...
    ./ max(max(abs(b(:)), [], 'omitnan'), eps);

errG = rel(bal.g, gCheck);
errF = rel(bal.f, fCheck);
errC = rel(bal.c, cCheck);
errGradAlpha = rel(bal.gradN_alpha, gradNAlphaCheck);
errGradAlphaSq = rel(bal.gradN_alpha_sq, gradNAlphaSqCheck);
errCross = rel(bal.cross_term, crossTermCheck);

fieldsToCheck = {'g','c','f','B','BoverBN','b_dot_gradN_theta', ...
    'gradN_alpha_sq','cross_term','BR','Bphi','BZ'};
allFinite = true;
fprintf('Miller ballooning coefficient diagnostic\n');
fprintf('  theta0                  : %.12g\n', bal.theta0);
fprintf('  nTheta                  : %d\n', numel(theta));
fprintf('  s_hat, alpha            : %.12g, %.12g\n', bal.s_hat, bal.alpha);
fprintf('  pprime                  : %.12e\n', bal.pprime);
fprintf('  FFprime                 : %.12e\n', bal.FFprime);
fprintf('  dpbar_dpsin             : %.12e\n', bal.dpbar_dpsin);

fprintf('\nFinite-value checks\n');
for k = 1:numel(fieldsToCheck)
    name = fieldsToCheck{k};
    val = bal.(name);
    ok = all(isfinite(val(:)));
    allFinite = allFinite && ok;
    fprintf('  %-22s : %d\n', name, ok);
end

fprintf('\nCoefficient self-consistency\n');
fprintf('  relerr g                : %.3e\n', errG);
fprintf('  relerr f                : %.3e\n', errF);
fprintf('  relerr c                : %.3e\n', errC);
fprintf('  relerr gradN_alpha      : %.3e\n', errGradAlpha);
fprintf('  relerr gradN_alpha_sq   : %.3e\n', errGradAlphaSq);
fprintf('  relerr cross_term       : %.3e\n', errCross);

fprintf('\nUseful ranges\n');
fprintf('  g min/max               : %.6e, %.6e\n', min(bal.g), max(bal.g));
fprintf('  f min/max               : %.6e, %.6e\n', min(bal.f), max(bal.f));
fprintf('  c min/max               : %.6e, %.6e\n', min(bal.c), max(bal.c));
fprintf('  B/B0 min/max            : %.6e, %.6e\n', min(bal.BoverBN), max(bal.BoverBN));
fprintf('  b dot gradN theta range : %.6e, %.6e\n', ...
    min(bal.b_dot_gradN_theta), max(bal.b_dot_gradN_theta));
fprintf('  gradN alpha sq range    : %.6e, %.6e\n', ...
    min(bal.gradN_alpha_sq), max(bal.gradN_alpha_sq));

periodicScalarFields = {'BoverBN','b_dot_gradN_theta','B','BR','Bphi','BZ'};
periodicVectorFields = {'bvec','gradN_psin','gradN_theta','gradN_zeta','gradN_H'};
nonperiodicScalarFields = {'g','f','c','gradN_alpha_sq','cross_term'};
nonperiodicVectorFields = {'gradN_alpha'};

fprintf('\nExpected periodic geometry/metric last-sample-to-first differences\n');
for k = 1:numel(periodicScalarFields)
    name = periodicScalarFields{k};
    jump = endpoint_jump(bal.(name));
    scale = max(max(abs(bal.(name)(:)), [], 'omitnan'), eps);
    fprintf('  %-22s : abs %.3e, rel %.3e\n', name, jump, jump/scale);
end
for k = 1:numel(periodicVectorFields)
    name = periodicVectorFields{k};
    jump = endpoint_jump(bal.(name));
    scale = max(vecnorm(bal.(name), 2, 2), [], 'omitnan');
    fprintf('  %-22s : abs %.3e, rel %.3e\n', name, jump, jump/max(scale, eps));
end

fprintf('\nAllowed non-periodic ballooning last-sample-to-first differences\n');
for k = 1:numel(nonperiodicScalarFields)
    name = nonperiodicScalarFields{k};
    jump = endpoint_jump(bal.(name));
    scale = max(max(abs(bal.(name)(:)), [], 'omitnan'), eps);
    fprintf('  %-22s : abs %.3e, rel %.3e\n', name, jump, jump/scale);
end
for k = 1:numel(nonperiodicVectorFields)
    name = nonperiodicVectorFields{k};
    jump = endpoint_jump(bal.(name));
    scale = max(vecnorm(bal.(name), 2, 2), [], 'omitnan');
    fprintf('  %-22s : abs %.3e, rel %.3e\n', name, jump, jump/max(scale, eps));
end

tol = 1e-10;
passed = allFinite ...
    && errG < tol ...
    && errF < tol ...
    && errC < tol ...
    && errGradAlpha < tol ...
    && errGradAlphaSq < tol ...
    && errCross < tol;

fprintf('\nOverall diagnostic pass     : %d\n', passed);

thetaClosed = [theta; 2*pi];

figure('Name', 'Miller periodic geometry and metric diagnostic');

subplot(2, 2, 1);
plot(thetaClosed, close_periodic(bal.BoverBN), 'LineWidth', 1.5);
hold on;
plot(thetaClosed, close_periodic(bal.b_dot_gradN_theta), 'LineWidth', 1.5);
grid on;
xlabel('\theta');
ylabel('periodic scalar');
legend('B/B_N', 'b\cdot\nabla_N\theta', 'Location', 'best');
title('Periodic scalar geometry');

subplot(2, 2, 2);
plot(thetaClosed, close_periodic(vecnorm(bal.gradN_psin, 2, 2)), 'LineWidth', 1.5);
hold on;
plot(thetaClosed, close_periodic(vecnorm(bal.gradN_theta, 2, 2)), 'LineWidth', 1.5);
plot(thetaClosed, close_periodic(vecnorm(bal.gradN_zeta, 2, 2)), 'LineWidth', 1.5);
plot(thetaClosed, close_periodic(vecnorm(bal.gradN_H, 2, 2)), 'LineWidth', 1.5);
grid on;
xlabel('\theta');
legend('|\nabla_N\psi_N|', '|\nabla_N\theta|', ...
    '|\nabla_N\zeta|', '|\nabla_N H|', 'Location', 'best');
title('Periodic vector-field norms');

subplot(2, 2, 3);
plot(thetaClosed, close_periodic(bal.BR), 'LineWidth', 1.5);
hold on;
plot(thetaClosed, close_periodic(bal.Bphi), 'LineWidth', 1.5);
plot(thetaClosed, close_periodic(bal.BZ), 'LineWidth', 1.5);
grid on;
xlabel('\theta');
legend('B_R', 'B_\phi', 'B_Z', 'Location', 'best');
title('Periodic magnetic-field components');

subplot(2, 2, 4);
plot([bal.R; bal.R(1)], [bal.Z; bal.Z(1)], 'LineWidth', 1.8);
axis equal;
grid on;
xlabel('R');
ylabel('Z');
title('Miller surface on PEST \theta');

figure('Name', 'Miller non-periodic ballooning coefficient diagnostic');

subplot(2, 2, 1);
plot(theta, bal.g, 'LineWidth', 1.5);
hold on;
plot(theta, bal.f, 'LineWidth', 1.5);
plot(theta, bal.c, 'LineWidth', 1.5);
grid on;
xlabel('\theta');
ylabel('coefficient');
legend('g', 'f', 'c', 'Location', 'best');
title('Non-periodic ballooning coefficient');

subplot(2, 2, 2);
plot(theta, bal.gradN_alpha(:,1), 'LineWidth', 1.5);
hold on;
plot(theta, bal.gradN_alpha(:,2), 'LineWidth', 1.5);
plot(theta, bal.gradN_alpha(:,3), 'LineWidth', 1.5);
grid on;
xlabel('\theta');
legend('R', '\phi', 'Z', 'Location', 'best');
title('\nabla_N\alpha components');

subplot(2, 2, 3);
plot(theta, bal.gradN_alpha_sq, 'LineWidth', 1.5);
hold on;
plot(theta, bal.cross_term, 'LineWidth', 1.5);
grid on;
xlabel('\theta');
legend('|\nabla_N\alpha|^2', '(b\times\nabla_N H)\cdot\nabla_N\alpha', ...
    'Location', 'best');
title('Non-periodic drive-related quantities');

subplot(2, 2, 4);
plot(theta, vecnorm(bal.gradN_alpha, 2, 2), 'LineWidth', 1.5);
grid on;
xlabel('\theta');
ylabel('|\nabla_N\alpha|');
title('Shear-driven secular growth');

function yClosed = close_periodic(y)
    yClosed = [y; y(1,:)];
end

function jump = endpoint_jump(y)
    jump = norm(y(end,:) - y(1,:));
end