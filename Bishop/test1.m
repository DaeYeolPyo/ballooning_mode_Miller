clc
clear
close all
%% Read GEQDSK and select flux surfaces
eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
%eq = read_geqdsk('./geqdsk_PT0.6');

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40); hold on;

psins = [0.2 0.4 0.6 0.8 1.0];
for k = 1:numel(psins)
    s = extract_flux_surface(eq, psins(k));
    plot(s.R, s.Z, 'LineWidth', 2);
end

plot(eq.rmaxis, eq.zmaxis, 'ko', 'MarkerFaceColor', 'y');
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Selected flux surfaces');

%% Extract a single flux surface and calculate tangent/normal vectors
surf = extract_flux_surface(eq, 0.5, ...
    'UniformTheta', true, ...
    'NTheta', 256, ...
    'ThetaCenter', 'axis');

geom = compute_tangent_normal(surf.R, surf.Z, ...
    'Closed', true, ...
    'NormalDirection', 'outward', ...
    'Method', 'Central');

figure;
plot(geom.R, geom.Z, 'k-', 'LineWidth', 1.5); hold on;
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Tangent and normal vectors on flux surface');

quiver(geom.R, geom.Z, geom.tR, geom.tZ, 0.15, ...
    'b', 'LineWidth', 1.2);
quiver(geom.R, geom.Z, geom.nR, geom.nZ, 0.15, ...
    'r', 'LineWidth', 1.2);

legend('Flux surface', 'Tangent', 'Normal');
grid on;

%% Calculate curvature vector
curv = compute_curvature(surf.R, surf.Z, ...
    'Closed', true, ...
    'NormalDirection', 'outward', ...
    'Method', 'central');

figure;
plot(curv.l, curv.kappa, 'LineWidth', 1.5);
xlabel('l [m]');
ylabel('\kappa [1/m]');
title('Curvature along flux surface');
grid on;

figure;
plot(curv.l, curv.Rc, 'LineWidth', 1.5);
xlabel('l [m]');
ylabel('R_c [m]');
title('Curvature radius along flux surface');
grid on;

%% Construct local (l, rho) coordinate around the flux surface
rho = linspace(-0.005, 0.005, 50);
loc = build_local_rho_coordinate(surf.R, surf.Z, rho, ...
    'Closed', true, ...
    'NormalDirection', 'outward');

figure; hold on;
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Local (l,\rho) coordinates around a flux surface');

% reference surface
plot(loc.R0, loc.Z0, 'k-', 'LineWidth', 2);

% rho = const lines
for j = 1:numel(loc.rho)
    plot(loc.R(:,j), loc.Z(:,j), 'b-');
end

% sample l = const connectors
idx = round(linspace(1, numel(loc.l), 50));
for i = idx
    plot(loc.R(i,:), loc.Z(i,:), 'r-');
end

%% Evaluate u(l)
uang = compute_u_angle(surf.R, surf.Z, ...
    'Closed', true, ...
    'Method', 'gradient', ...
    'NormalDirection', 'outward', ...
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

%% Evaluate h0(l)
h = compute_h0(surf.R, surf.Z, ...
    'Closed', true, ...
    'Method', 'gradient');

figure;
plot(h.l, h.h0_from_R, 'LineWidth', 1.5, 'DisplayName', 'h_0 = R/R_0'); hold on;
plot(h.l, h.h0_from_integral, '--', 'LineWidth', 1.5, ...
    'DisplayName', 'h_0 = 1 + (1/R_0)\int cos(u) dl');
xlabel('l [m]');
ylabel('h_0');
title('Check of h_0(l)');
legend show;
grid on;