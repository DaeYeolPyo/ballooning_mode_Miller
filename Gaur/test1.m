clc
close all
clear

%eq = read_geqdsk('../Bishop/geqdsk_PT0.6');
eq = read_geqdsk('../Bishop/curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');

%----------------------------------
% 설정
%----------------------------------
psin_list   = linspace(0.05, 0.95, 30);
theta_lines = linspace(0, 2*pi, 300);
theta_lines(end) = [];   % 0과 2pi 중복 제거

JacMode = 1;             % 1 = equal-arc-length, 4 = PEST
NTheta  = 256;

npsi = length(psin_list);
ntheta_line = length(theta_lines);

R_all = cell(npsi,1);
Z_all = cell(npsi,1);
theta_all = cell(npsi,1);

%----------------------------------
% 1. 각 flux surface에서 straight-field-line theta 계산
%----------------------------------
for i = 1:npsi
    psin = psin_list(i);

    cont = eq_straight_fieldline_theta(eq, psin, ...
        'JacMode', JacMode, ...
        'NTheta', NTheta);

    R_all{i} = cont.R(:);
    Z_all{i} = cont.Z(:);
    theta_all{i} = cont.theta_coord(:);
end

%----------------------------------
% 2. flux surfaces 그리기
%----------------------------------
figure; hold on;

for i = 1:npsi
    plot(R_all{i}, Z_all{i}, 'k-', 'LineWidth', 0.5);
end

%----------------------------------
% 3. equal-theta lines 그리기
%----------------------------------
for it = 1:ntheta_line
    theta0 = theta_lines(it);

    R_line = nan(npsi,1);
    Z_line = nan(npsi,1);

    for i = 1:npsi
        th = theta_all{i};
        R  = R_all{i};
        Z  = Z_all{i};

        % periodic extension
        th_ext = [th; th(1) + 2*pi];
        R_ext  = [R; R(1)];
        Z_ext  = [Z; Z(1)];

        R_line(i) = interp1(th_ext, R_ext, theta0, 'linear');
        Z_line(i) = interp1(th_ext, Z_ext, theta0, 'linear');
    end

    plot(R_line, Z_line, 'r-', 'LineWidth', 1.5);
end

%----------------------------------
% 4. 그림 설정
%----------------------------------
grid on;

xlabel('R');
ylabel('Z');

if JacMode == 4
    coordName = 'PEST';
else
    coordName = sprintf('JacMode = %d', JacMode);
end

title(['Equal-\theta lines in ', coordName, ' coordinate']);