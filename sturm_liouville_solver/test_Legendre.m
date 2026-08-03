clc
clear
close all

%% Numerical parameters
ngrid = 200;

nmode = 100;
numplot = 5;

%% Define g, c, f
% Legendre polynomial test problem
% g = sin(theta), c = 0, f = sin(theta)
% will give X = P_l(cos(theta))
% lambda = -l(l+1) (l = 0, 1, 2, ...)
thetaGrid = linspace(0, pi, ngrid);
g = sin(thetaGrid);
c = zeros(1, ngrid);
f = sin(thetaGrid);

%% Construct and solve eigenvalue problem
% Construct FEM matrices
[Gmat, Cmat, Fmat, thetaDof] = construct_matrix(thetaGrid, g, c, f);
% Solve eigenvalue
[lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat, 'natural');

%% Graphics
mode = 0:nmode-1;
lambdaExact = -mode.*(mode + 1);

figure(1);
plot(mode, lambda(1:numel(mode)), 'bo--', 'DisplayName', 'Numerical');
hold on;
plot(mode, lambdaExact, 'r-', 'DisplayName', 'Analytic');
hold off;
legend;

for i = 1:numplot
    ell = i - 1;

    Xnum = X(:, i);
    Xana = legendreP(ell, cos(thetaDof));
    Xana = Xana(:);

    % Normalization w.r.t Fmat such that X^T F X = 1
    Xnum = Xnum / sqrt(real(Xnum.' * Fmat * Xnum));
    Xana = Xana / sqrt(real(Xana.' * Fmat * Xana));

    % Sign synchronyzation
    if real(Xana.' * Fmat * Xnum) < 0
        Xnum = -Xnum;
    end

    figure;
    plot(thetaDof, Xnum, 'b--', 'DisplayName', 'Numerical');
    hold on;
    plot(thetaDof, Xana, 'r-', 'DisplayName', 'Analytic');
    hold off;
    legend;
end