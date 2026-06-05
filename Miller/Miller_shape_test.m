clc
clear
close all

%% Read GEQDSK and select flux surfaces
eq = read_geqdsk('./geqdsk_PT0.7');

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40); hold on;

psins = [0.3, 0.5, 0.9];
for k = 1:numel(psins)
    s = extract_flux_surface(eq, psins(k));
    plot(s.R, s.Z, 'LineWidth', 2);
end

plot(eq.rbbbs, eq.zbbbs, '-k', 'LineWidth', 3);

plot(eq.rmaxis, eq.zmaxis, 'ko', 'MarkerFaceColor', 'y');
xlabel('R [m]');
ylabel('Z [m]');
title('Selected flux surfaces');

%% Fit the flux surface at psiN = 0.77 [Miller(1998)]
psiN = 0.77;
param = fit_Miller(eq, psiN);

figure;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40); hold on;
s = extract_flux_surface(eq, psiN);
[R, Z, u, p] = DshapeParam(param, 300);
plot(s.R, s.Z, '-k');
plot(R, Z, '-r');
xlabel('R [m]');
ylabel('Z [m]');
title('Fitted Miller shape');