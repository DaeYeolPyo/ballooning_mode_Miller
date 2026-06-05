clc
clear
close all
%% Read GEQDSK and select flux surfaces
eq = read_geqdsk('../Bishop/curavture_analysis/geqdsk/scan_B2.5_C18_G08_H11.geqdsk');
%eq = read_geqdsk('./geqdsk_PT0.6');

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40); hold on;

psins = [0.3, 0.5, 0.9];
for k = 1:numel(psins)
    s = extract_flux_surface(eq, psins(k));
    R_Z0 = find_R_at_Z0_from_surface(s);
    plot(s.R, s.Z, 'LineWidth', 2);
    plot(R_Z0, zeros(size(R_Z0)), 'ro', 'MarkerFaceColor', 'r');
    yline(0, '--');
end

plot(eq.rmaxis, eq.zmaxis, 'ko', 'MarkerFaceColor', 'y');
axis equal;
xlabel('R [m]');
ylabel('Z [m]');
title('Selected flux surfaces');

%% Fit the flux surfaces into parametrized shape
ntheta = 300;

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40); hold on;

for k = 1:numel(psins)
    s = extract_flux_surface(eq, psins(k));
    coeffs = fit_Cshape(s);
    plot(s.R, s.Z, 'LineWidth', 2);
    [fitR, fitZ] = CshapeParam(coeffs, ntheta);
    plot(fitR, fitZ, 'LineStyle', '--');
end

%% Calculate derivative of shaping coefficients
psiN = 0.9;

coeffs = Cshape_metrics_bpfit(eq, psiN, ntheta);
[Rr, Zr] = r_derv(coeffs, ntheta);
[Rt, Zt] = theta_derv(coeffs, ntheta);

Bp = Bpol(eq, psiN, ntheta, coeffs);
figure;

subplot(1, 2, 1);
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40);
hold on;
surf = extract_flux_surface(eq, psiN);
plot(surf.R, surf.Z, 'LineWidth', 2);
[R, Z] = CshapeParam(coeffs, ntheta);
plot(R(1), Z(1), 'ro');
plot(R, Z, 'LineStyle', '--');
legend('', 'Flux surface', '', 'Fitted surface')

subplot(1, 2, 2);
Bp_EFIT = Bpol_from_EFIT(eq, R, Z);
arc = cumsum([0, hypot(diff(R), diff(Z))]);
arcTheta = 2*pi*arc/arc(end);

plot(arcTheta, Bp);
hold on;
%plot(arcTheta, Bp_EFIT);
xlabel('Equal arc-length angle');
ylabel('B_p');
%legend('C-shape local model', 'EFIT', 'Location', 'best');
