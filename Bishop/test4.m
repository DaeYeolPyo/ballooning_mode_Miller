clc
clear
close all
%% Read GEQDSK and extract a single flux surface
%eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
eq = read_geqdsk('./geqdsk_PT0.6');
psiN = 0.9;
NTheta = 1024;

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40);
hold on

% Bishop report convention:
%   l increases clockwise.
%   Bishop normal is n = (sin u, -cos u).
% JacMode = 1 builds an equal-arc-length poloidal coordinate, then the
% contour is resampled on uniform l for cleaner derivatives.
cont = eq_straight_fieldline_theta(eq, psiN, ...
    'JacMode', 1, 'NTheta', NTheta, ...
    'Direction', 'clockwise');

surf_l = resample_closed_curve_arclength(cont.R, cont.Z, NTheta, ...
    'Direction', 'clockwise');

plot([surf_l.R; surf_l.R(1)], ...
    [surf_l.Z; surf_l.Z(1)], ...
    '-k', 'LineWidth', 2);
plot(eq.rmaxis, eq.zmaxis, 'ko', 'MarkerFaceColor', 'y');
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Selected flux surfaces');

%% Check
R = surf_l.R;
Z = surf_l.Z;

geom = compute_tangent_normal(R, Z, ...
    'Closed', true, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'Method', 'Central');

figure;
plot(geom.R, geom.Z, 'k-', 'LineWidth', 1.5); hold on;
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Clockwise tangent and Bishop normal on flux surface');

quiver(geom.R, geom.Z, geom.tR, geom.tZ, 0.15, ...
    'b', 'LineWidth', 1.2);
quiver(geom.R, geom.Z, geom.nR, geom.nZ, 0.15, ...
    'r', 'LineWidth', 1.2);

legend('Flux surface', 'Tangent', 'Bishop normal');
grid on;

%% Calculate curvature vector
curv = compute_curvature(R, Z, ...
    'Closed', true, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'Method', 'central');

figure;
plot(curv.l, curv.kappa_signed, 'LineWidth', 1.5); hold on;
plot(curv.l, curv.kappa, '--', 'LineWidth', 1.2);
xlabel('l [m]');
ylabel('\kappa [1/m]');
title('Curvature along flux surface');
legend('\kappa signed (Bishop normal)', '|\kappa|');
grid on;

figure;
plot(curv.l, curv.Rc_signed, 'LineWidth', 1.5); hold on;
plot(curv.l, curv.Rc, '--', 'LineWidth', 1.2);
xlabel('l [m]');
ylabel('R_c [m]');
title('Curvature radius along flux surface');
legend('R_c signed', '|R_c|');
grid on;

%% Evaluate u(l)
uang = compute_u_angle(R, Z, ...
    'Closed', true, ...
    'Method', 'central', ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'Unwrap', true);

figure;
plot(uang.l, uang.u_raw, 'o-', 'DisplayName', 'u raw'); hold on;
plot(uang.l, uang.u, '-', 'LineWidth', 1.5, 'DisplayName', 'u unwrapped');
xlabel('l [m]');
ylabel('u [rad]');
title('Tangent angle u(l)');
legend show;
grid on;

figure;
subplot(2,1,1);
plot(uang.l, uang.tR, 'LineWidth', 1.5); hold on;
plot(uang.l, cos(uang.u), '--', 'LineWidth', 1.5);
ylabel('dR/dl , cos(u)');
legend('t_R','cos(u)');
grid on;

subplot(2,1,2);
plot(uang.l, uang.tZ, 'LineWidth', 1.5); hold on;
plot(uang.l, sin(uang.u), '--', 'LineWidth', 1.5);
ylabel('dZ/dl , sin(u)');
xlabel('l [m]');
legend('t_Z','sin(u)');
grid on;

%% Extend one-period geometry for ballooning representation
base = geom;
base.u = uang.u;
base.theta_ballooning = 2*pi * base.l / base.L;

balloon = extend_ballooning_theta(base, ...
    'NPeriodsEachSide', 3, ...
    'Direction', 'clockwise');

figure;
plot(balloon.theta, balloon.R, 'LineWidth', 1.2); hold on;
plot(balloon.theta, balloon.Z, 'LineWidth', 1.2);
xlabel('\theta');
ylabel('R, Z [m]');
title('Finite window of ballooning theta extension');
legend('R(\theta)', 'Z(\theta)');
grid on;

fprintf('Bishop geometry check: direction=%s, normal=%s, L=%.8g, u span=%.8g rad\n', ...
    geom.direction, geom.normalDirection, geom.L, uang.u(end)-uang.u(1));
fprintf('Ballooning theta window: [%.3f, %.3f], points=%d\n', ...
    min(balloon.theta), max(balloon.theta), numel(balloon.theta));
