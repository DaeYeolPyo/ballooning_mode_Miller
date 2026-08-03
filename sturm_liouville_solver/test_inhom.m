clc
clear
close all

%% Numerical parameters
ngrid = 200;

c0 = 0;

nmode = 100;
numplot = 5;

%% Define g, c, f
% Inhomogeneous test problem
% g = 1 + theta, c = c0/(1 + theta), f = 1/(1 + theta)
% will give X = sin(n*pi*ln(1 + theta)/ln(2))
% lambda = c0 - (n*pi/ln(2))^2 (n = 1, 2, ...)
thetaGrid = linspace(0, 1, ngrid);
g = 1 + thetaGrid;
c = c0./(1 + thetaGrid);
f = 1./(1 + thetaGrid);

%% Construct and solve eigenvalue problem
% Construct FEM matrices
[Gmat, Cmat, Fmat, thetaDof] = construct_matrix(thetaGrid, g, c, f);
% Solve eigenvalue
[lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat);

%% Graphics
mode = 1:nmode;
figure(1);
plot(mode, lambda(1:numel(mode)), 'bo--', 'DisplayName', 'Numerical');
hold on;
plot(mode, c0 - (mode*pi/log(2)).^2, 'r-', 'DisplayName', 'Analytic');
hold off;
legend;

for i = 1:numplot
    ell = i - 1;

    Xnum = X(:, i);
    Xana = sin(i*pi*log(1 + thetaDof)/log(2));
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