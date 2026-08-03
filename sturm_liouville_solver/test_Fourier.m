clc
clear
close all

%% Numerical parameters
ngrid = 200;

c0 = 3;

nmode = 30;
numplot = 5;

%% Define g, c, f
% General Fourier harmonics problem
% g = 1, c = c0, f = 1
% will give X = exp(im*theta)
% lambda = c0 - m^2 (m = 0, 1, 2, ...)
thetaGrid = linspace(0, 2*pi, ngrid);
g = ones(1, ngrid);
c = c0.*ones(1, ngrid);
f = ones(1, ngrid);

%% Construct and solve eigenvalue problem
% Construct FEM matrices
[Gmat, Cmat, Fmat, thetaDof] = construct_matrix(thetaGrid, g, c, f);
% Solve eigenvalue
[lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat, 'periodic');

%% Graphics
mode = 1:nmode;
m = ceil((0:nmode-1)/2);
lambdaExact = c0 - m.^2;

figure(1);
plot(mode, lambda(1:nmode), 'bo--', ...
    'DisplayName', 'Numerical');
hold on;
plot(mode, lambdaExact, 'r-', ...
    'DisplayName', 'Analytic');
hold off;
legend;

for i = 1:numplot
    mi = m(i);
    Xnum = X(:, i);

    if mi == 0
        % m = 0: 상수 모드
        Xana = ones(size(thetaDof));
    else
        % 수치해는 cos(mi*theta), sin(mi*theta)의
        % 임의 선형결합으로 나올 수 있다.
        Phi = [cos(mi*thetaDof), sin(mi*thetaDof)];

        % 수치해에 가장 가까운 Fourier 해석해의 위상/진폭을 찾음
        coeff = (Phi.' * Fmat * Phi) \ (Phi.' * Fmat * Xnum);
        Xana = Phi * coeff;
    end

    % Fmat 기준 정규화
    Xnum = Xnum / sqrt(real(Xnum.' * Fmat * Xnum));
    Xana = Xana / sqrt(real(Xana.' * Fmat * Xana));

    % 부호 통일
    if real(Xana.' * Fmat * Xnum) < 0
        Xnum = -Xnum;
    end

    figure;
    plot(thetaDof, Xnum, 'b--', 'DisplayName', 'Numerical');
    hold on;
    plot(thetaDof, Xana, 'r-', 'DisplayName', 'Analytic');
    hold off;
    title(sprintf('Fourier mode %d, m = %d', i, mi));
    legend;
end