clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, '..', 'Gaur'));

fprintf('standard s-alpha solver check\n');

thetaB = 5*pi;
nGrid = 401;
cases = [
    0.0, 0.0
    1.0, 1.0
    2.0, 2.0
    3.0, 4.0
    5.0, 7.0
];

fprintf('  s_hat   alpha      solver lambda      direct lambda       abs diff\n');
for i = 1:size(cases, 1)
    shat = cases(i, 1);
    alpha = cases(i, 2);
    bal = standard_salpha_bal(shat, alpha, 257);

    sol = solve_ballooning_eigenvalue(bal, ...
        'ThetaB', thetaB, ...
        'N', nGrid, ...
        'NEigs', 6, ...
        'Plot', false);

    lambdaDirect = direct_standard_salpha_lambda(shat, alpha, thetaB, nGrid);

    fprintf('  %5.2f  %6.2f  % .12e  % .12e  %.3e\n', ...
        shat, alpha, real(sol.lambda), lambdaDirect, ...
        abs(real(sol.lambda) - lambdaDirect));
end

fprintf('\nconstant-coefficient analytic limit\n');
bal = standard_salpha_bal(0.0, 0.0, 257);
sol = solve_ballooning_eigenvalue(bal, ...
    'ThetaB', thetaB, ...
    'N', nGrid, ...
    'NEigs', 6, ...
    'Plot', false);
lambdaExact = -(pi/(2*thetaB))^2;
fprintf('  exact  lambda : %.12e\n', lambdaExact);
fprintf('  solver lambda : %.12e\n', real(sol.lambda));
fprintf('  abs err       : %.3e\n', abs(real(sol.lambda) - lambdaExact));

function bal = standard_salpha_bal(shat, alpha, ntheta)
    theta = linspace(0, 2*pi, ntheta+1).';
    theta(end) = [];

    bal = struct();
    bal.theta = theta;
    bal.q = 0.0;
    bal.dq_dpsin = shat;
    bal.dpbar_dpsin = 1.0;
    bal.BoverBN = ones(ntheta, 1);
    bal.b_dot_gradN_theta = ones(ntheta, 1);
    bal.bvec = repmat([0, 0, 1], ntheta, 1);
    bal.gradN_psin = repmat([1, 0, 0], ntheta, 1);
    bal.gradN_theta = zeros(ntheta, 3);
    bal.gradN_zeta = repmat([0, 1, 0], ntheta, 1);
    bal.gradN_H = [0.5*alpha*cos(theta), 0.5*alpha*sin(theta), zeros(ntheta, 1)];
    bal.local_shear_model = 'miller-s-alpha';
    bal.local_shear_periodic_alpha = alpha*ones(ntheta, 1);
    bal.theta0 = 0.0;
end

function lambda = direct_standard_salpha_lambda(shat, alpha, thetaB, nGrid)
    nGrid = round(nGrid);
    if mod(nGrid, 2) == 0
        nGrid = nGrid + 1;
    end

    theta = linspace(-thetaB, thetaB, nGrid).';
    dtheta = theta(2) - theta(1);

    thetaHalf = 0.5*(theta(1:end-1) + theta(2:end));
    h = shat*theta - alpha*sin(theta);
    hHalf = shat*thetaHalf - alpha*sin(thetaHalf);

    g = 1 + h.^2;
    f = g;
    c = alpha*(cos(theta) + h.*sin(theta));
    gHalf = 1 + hHalf.^2;

    nInterior = nGrid - 2;
    lowerK = zeros(nInterior, 1);
    diagK = zeros(nInterior, 1);
    upperK = zeros(nInterior, 1);

    for m = 1:nInterior
        gMinus = gHalf(m);
        gPlus = gHalf(m+1);
        diagK(m) = c(m+1) - (gMinus + gPlus)/dtheta^2;
        if m > 1
            lowerK(m) = gMinus/dtheta^2;
        end
        if m < nInterior
            upperK(m) = gPlus/dtheta^2;
        end
    end

    K = diag(diagK) + diag(upperK(1:end-1), 1) + diag(lowerK(2:end), -1);
    M = diag(f(2:end-1));
    lambda = max(real(eig(K, M)));
end
