clc
clear
close all

%% Numerical parameters
ngrid = 200;
m = 1; % Order
R = 1; % Radius

nmode = 100;
numplot = 5;

%% Define g, c, f
% Bessel function test problem
% g = theta^(2m+1), c = 0, f = theta^(2m+1)
% will give X = J_m(j_mn theta/R)
% lambda = -(j_mn/R)^2
thetaGrid = linspace(0, R, ngrid);
g = thetaGrid.^(2*m+1);
c = zeros(size(thetaGrid));
f = thetaGrid.^(2*m+1);

%% Construct and solve eigenvalue problem
% Construct FEM matrices
[Gmat, Cmat, Fmat, thetaDof] = construct_matrix(thetaGrid, g, c, f);
% Solve eigenvalue
[lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat, 'natural-dirichlet');

%% Graphics
mode = 1:nmode;
jmn = zeros(size(mode));

zScan = linspace(1e-8, (nmode + m + 5)*pi, 10000);
values = besselj(m, zScan);
crossings = find(values(1:end-1).*values(2:end) < 0);

for i = 1:nmode
    bracket = zScan(crossings(i):crossings(i)+1);
    jmn(i) = fzero(@(z) besselj(m, z), bracket);
end

lambdaExact = -(jmn/R).^2;

figure(1);
plot(mode, lambda(1:numel(mode)), 'bo--', 'DisplayName', 'Numerical');
hold on;
plot(mode, lambdaExact, 'r-', 'DisplayName', 'Analytic');
hold off;
legend;

for i = 1:numplot
    ell = mode(i);

    Xnum = thetaDof.^m.*X(:, i);
    Xana = besselj(m, jmn(i)*thetaDof/R);
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