clc
clear
close all

%% Numerical parameters
ngrid = 200;
m = 1; % Order

nmode = 100;
numplot = 5;

%% Define g, c, f
% Associated Legendre function test problem
% g = sin^(2m+1)(theta), c = -m(m+1)sin^(2m+1)(theta), f = sin^(2m+1)(theta)
% will give X = P_l^m(cos(theta)).sin^(m)(theta)
% lambda = -l(l+1) (l >= |m|)
thetaGrid = linspace(0, pi, ngrid);
sinArr = sin(thetaGrid).^(2*m+1);
g = sinArr;
c = -m*(m+1).*sinArr;
f = sinArr;

%% Construct and solve eigenvalue problem
% Construct FEM matrices
[Gmat, Cmat, Fmat, thetaDof] = construct_matrix(thetaGrid, g, c, f);
% Solve eigenvalue
[lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat, 'natural');

%% Graphics
mode = m + (0:nmode-1);
lambdaExact = -mode.*(mode + 1);

figure(1);
plot(mode, lambda(1:numel(mode)), 'bo--', 'DisplayName', 'Numerical');
hold on;
plot(mode, lambdaExact, 'r-', 'DisplayName', 'Analytic');
hold off;
xlabel('l')
ylabel('\lambda=-l(l+1)')
legend;

for i = 1:numplot
    ell = mode(i);

    Xnum = X(:, i).*(sin(thetaDof).^m);
    Pall = legendre(ell, cos(thetaDof));
    Xana = Pall(m + 1, :).';
    Xana = Xana(:);

    % Normalization w.r.t Fmat such that X^T F X = 1
    Xnum = Xnum / sqrt(real(Xnum.' * Fmat * Xnum));
    Xana = Xana / sqrt(real(Xana.' * Fmat * Xana));

    % Sign synchronyzation
    if real(Xana.' * Fmat * Xnum) < 0
        Xnum = -Xnum;
    end

%     figure;
    plot(thetaDof, Xnum, '--', 'DisplayName', sprintf('P_%d^1(cos\\theta)', ell));
    hold on;
%     plot(thetaDof, Xana, 'r-', 'DisplayName', 'Analytic');
%     hold off;
    legend;
end
xlabel('\theta')