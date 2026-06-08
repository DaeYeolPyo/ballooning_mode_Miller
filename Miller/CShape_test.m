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

%% Calculate derivative of shaping coefficients and check poloidal field
psiN = 0.9;

coeffs_geom = Cshape_metrics(eq, psiN);
coeffs_bpfit = Cshape_metrics_bpfit(eq, psiN, ntheta);
[Rr, Zr] = r_derv(coeffs_bpfit, ntheta);
[Rt, Zt] = theta_derv(coeffs_bpfit, ntheta);

Bp_geom = Bpol(eq, psiN, ntheta, coeffs_geom);
Bp_bpfit = Bpol(eq, psiN, ntheta, coeffs_bpfit);
figure;

subplot(1, 2, 1);
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40);
hold on;
surf = extract_flux_surface(eq, psiN);
plot(surf.R, surf.Z, 'LineWidth', 2);
[R, Z] = CshapeParam(coeffs_bpfit(:, 1), ntheta);
plot(R(1), Z(1), 'ro');
plot(R, Z, 'LineStyle', '--');
legend('', 'Flux surface', '', 'Fitted surface')

subplot(1, 2, 2);
Bp_EFIT = Bpol_from_EFIT(eq, R, Z);
arc = cumsum([0, hypot(diff(R), diff(Z))]);
arcTheta = 2*pi*arc/arc(end);

plot(arcTheta, Bp_geom, '-k', 'LineWidth', 1.5);
hold on;
%plot(arcTheta, Bp_bpfit, '-b', 'LineWidth', 1.5);
plot(arcTheta, Bp_EFIT, '-r', 'LineWidth', 1.5);
grid on;
xlabel('\theta');
ylabel('B_p [T]');
legend('From fitting', 'From geqdsk', 'Location', 'best');

relerr = @(a,b) max(abs(a(:) - b(:)), [], 'omitnan') ...
    ./ max(max(abs(b(:)), [], 'omitnan'), eps);
fprintf('C-shape Bpol diagnostic at psiN = %.6g\n', psiN);
fprintf('  geometry-fit relerr max : %.6e\n', relerr(Bp_geom, Bp_EFIT));
fprintf('  Bp-fit relerr max       : %.6e\n', relerr(Bp_bpfit, Bp_EFIT));
fprintf('  EFIT Bp min/max         : %.6e, %.6e\n', ...
    min(Bp_EFIT, [], 'omitnan'), max(Bp_EFIT, [], 'omitnan'));
fprintf('  Bp-fit min/max          : %.6e, %.6e\n', ...
    min(Bp_bpfit, [], 'omitnan'), max(Bp_bpfit, [], 'omitnan'));
