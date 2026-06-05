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