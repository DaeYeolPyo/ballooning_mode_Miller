clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'salpha'));
addpath(fullfile(thisDir, '..', 'Gaur'));

fprintf('solve_ballooning_eigenvalue diagnostic\n');

%% Case 1: constant coefficients with known analytic eigenvalue.
% d2X/dtheta^2 + c X = lambda f X, X(-L)=X(L)=0.
% For g=1, c=0, f=1, the largest eigenvalue is
% lambda_1 = -(pi/(2L))^2.
L = 5*pi;
N = 401;
thetaBase = linspace(0, 2*pi, 129).';
thetaBase(end) = [];

balConst = struct();
balConst.theta = thetaBase;
balConst.g = ones(size(thetaBase));
balConst.c = zeros(size(thetaBase));
balConst.f = ones(size(thetaBase));
balConst.theta0 = 0;

solConst = solve_ballooning_eigenvalue(balConst, ...
    'ThetaB', L, ...
    'N', N, ...
    'NEigs', 6, ...
    'UseSparse', true, ...
    'Plot', false);

lambdaExact = -(pi/(2*L))^2;
lambdaConst = real(solConst.lambda);
lambdaConstRefined = real(solConst.lambda_refined);

fprintf('\nConstant-coefficient analytic check\n');
fprintf('  exact lambda             : %.12e\n', lambdaExact);
fprintf('  solver lambda            : %.12e\n', lambdaConst);
fprintf('  refined lambda           : %.12e\n', lambdaConstRefined);
fprintf('  abs err solver           : %.3e\n', abs(lambdaConst - lambdaExact));
fprintf('  abs err refined          : %.3e\n', abs(lambdaConstRefined - lambdaExact));
fprintf('  boundary X endpoints     : %.3e, %.3e\n', solConst.X(1), solConst.X(end));
fprintf('  used full eig            : %d\n', solConst.used_full_eig);

%% Case 2: same constant problem with full eig vs sparse eig consistency.
solConstFull = solve_ballooning_eigenvalue(balConst, ...
    'ThetaB', L, ...
    'N', N, ...
    'NEigs', 6, ...
    'UseSparse', false, ...
    'FullEigThreshold', inf, ...
    'Plot', false);

fprintf('\nSparse/full eig consistency check\n');
fprintf('  sparse lambda            : %.12e\n', real(solConst.lambda));
fprintf('  full lambda              : %.12e\n', real(solConstFull.lambda));
fprintf('  abs difference           : %.3e\n', abs(real(solConst.lambda) - real(solConstFull.lambda)));

%% Case 3: realistic Miller coefficient solve.
p = DshapeMillerParams();
p.delta = 0.7;
theta0 = 0.0;

balMiller = miller_ballooning_coefficients(p, ...
    'Theta0', theta0, ...
    'NTheta', 257, ...
    'NGeom', 1201);

solMiller = solve_ballooning_eigenvalue(balMiller, ...
    'Theta0', theta0, ...
    'ThetaB', L, ...
    'N', N, ...
    'NEigs', 6, ...
    'Plot', false);

residual = solMiller.K*solMiller.X(2:end-1) ...
    - solMiller.lambda*solMiller.M*solMiller.X(2:end-1);
residualNorm = norm(residual) ...
    / max(norm(solMiller.K*solMiller.X(2:end-1)), eps);

fprintf('\nMiller-coefficient solve check\n');
fprintf('  theta0                  : %.12g\n', theta0);
fprintf('  lambda                  : %.12e\n', real(solMiller.lambda));
fprintf('  lambda refined          : %.12e\n', real(solMiller.lambda_refined));
fprintf('  lambda-refined diff     : %.3e\n', abs(real(solMiller.lambda) - real(solMiller.lambda_refined)));
fprintf('  relative residual        : %.3e\n', residualNorm);
fprintf('  boundary X endpoints     : %.3e, %.3e\n', solMiller.X(1), solMiller.X(end));
fprintf('  g min/max                : %.6e, %.6e\n', min(solMiller.g), max(solMiller.g));
fprintf('  f min/max                : %.6e, %.6e\n', min(solMiller.f), max(solMiller.f));
fprintf('  c min/max                : %.6e, %.6e\n', min(solMiller.c), max(solMiller.c));
fprintf('  used full eig            : %d\n', solMiller.used_full_eig);

figure('Name', 'solve_ballooning_eigenvalue constant-coefficient check');

Xflip = interp1(solMiller.theta_b, solMiller.X, -solMiller.theta_b, 'linear', NaN);
validSym = isfinite(Xflip);
evenError = norm(solMiller.X(validSym) - Xflip(validSym)) ...
    / max(norm(solMiller.X(validSym)), eps);
oddError = norm(solMiller.X(validSym) + Xflip(validSym)) ...
    / max(norm(solMiller.X(validSym)), eps);

gFlip = interp1(solMiller.theta_b, solMiller.g, -solMiller.theta_b, 'linear', NaN);
fFlip = interp1(solMiller.theta_b, solMiller.f, -solMiller.theta_b, 'linear', NaN);
cFlip = interp1(solMiller.theta_b, solMiller.c, -solMiller.theta_b, 'linear', NaN);
validCoefSym = isfinite(gFlip) & isfinite(fFlip) & isfinite(cFlip);
gEvenError = norm(solMiller.g(validCoefSym) - gFlip(validCoefSym)) ...
    / max(norm(solMiller.g(validCoefSym)), eps);
fEvenError = norm(solMiller.f(validCoefSym) - fFlip(validCoefSym)) ...
    / max(norm(solMiller.f(validCoefSym)), eps);
cEvenError = norm(solMiller.c(validCoefSym) - cFlip(validCoefSym)) ...
    / max(norm(solMiller.c(validCoefSym)), eps);

fprintf('\nMiller symmetry check about theta_b = 0\n');
fprintf('  X even relative error     : %.3e\n', evenError);
fprintf('  X odd relative error      : %.3e\n', oddError);
fprintf('  g even relative error     : %.3e\n', gEvenError);
fprintf('  f even relative error     : %.3e\n', fEvenError);
fprintf('  c even relative error     : %.3e\n', cEvenError);


subplot(2, 2, 1);
plot(solConst.theta_b, solConst.X, 'LineWidth', 1.6);
grid on;
xlabel('\theta');
ylabel('X');
title('Constant-coefficient eigenfunction');

subplot(2, 2, 2);
plot(solConst.theta_b, solConst.g, 'LineWidth', 1.5);
hold on;
plot(solConst.theta_b, solConst.f, 'LineWidth', 1.5);
plot(solConst.theta_b, solConst.c, 'LineWidth', 1.5);
grid on;
xlabel('\theta');
legend('g', 'f', 'c', 'Location', 'best');
title('Constant coefficients on ballooning grid');

subplot(2, 2, 3);
plot(solConst.theta_b, solConst.X, 'LineWidth', 1.6);
hold on;
plot(solConst.theta_b, cos(pi*solConst.theta_b/(2*L)), '--', 'LineWidth', 1.3);
grid on;
xlabel('\theta');
legend('solver', 'analytic shape', 'Location', 'best');
title('Eigenfunction shape check');

subplot(2, 2, 4);
bar([lambdaExact, lambdaConst, lambdaConstRefined]);
grid on;
set(gca, 'XTickLabel', {'exact', 'solver', 'refined'});
ylabel('\lambda');
title('Eigenvalue comparison');

figure('Name', 'solve_ballooning_eigenvalue Miller-coefficient check');

subplot(2, 2, 1);
plot(solMiller.theta_b, solMiller.X, 'LineWidth', 1.6);
hold on;
plot(solMiller.theta_b, Xflip, '--', 'LineWidth', 1.2);
grid on;
xlabel('\theta');
ylabel('X');
legend('X(\theta)', 'X(-\theta)', 'Location', 'best');
title('Eigenfunction X');

subplot(2, 2, 2);
plot(solMiller.theta_b, solMiller.g, 'LineWidth', 1.4);
hold on;
plot(solMiller.theta_b, solMiller.f, 'LineWidth', 1.4);
plot(solMiller.theta_b, solMiller.c, 'LineWidth', 1.4);
grid on;
xlabel('\theta');
legend('g', 'f', 'c', 'Location', 'best');
title('Extended ballooning coefficients');

subplot(2, 2, 3);
semilogy(abs(residual), 'LineWidth', 1.4);
grid on;
xlabel('interior grid index');
ylabel('|AX - \lambda X|');
title('Discrete residual');