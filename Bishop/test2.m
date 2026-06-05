clc
clear
close all

eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H09.geqdsk');

geom = build_bishop_geometry(eq, 0.7);

figure;
plot(geom.l, geom.h0)
title('h0 = R/R0')

figure;
dudl = gradient(geom.u, geom.l);

plot(geom.l, abs(dudl), geom.l, geom.kappa)
legend('|du/dl|','kappa')

figure;
plot(geom.l, geom.Bp)
title('Poloidal field B_p')

figure;
subplot(3,1,1)
plot(geom.l, geom.I, 'LineWidth', 1.5)
ylabel('I')

subplot(3,1,2)
plot(geom.l, geom.Iprime, 'LineWidth', 1.5)
ylabel('I''')

subplot(3,1,3)
plot(geom.l, geom.pprime, 'LineWidth', 1.5)
ylabel('p''')
xlabel('l')
grid on

eq29 = compute_eq29(geom, 'L0Index', 1);
eq31 = compute_eq31(geom, eq29, 'UseMu0PprimeDrive', true);

figure;
subplot(3,1,1)
plot(eq29.l, eq29.integrand, 'LineWidth', 1.5)
ylabel('integrand')
grid on

subplot(3,1,2)
plot(eq29.l, eq29.integral_term, 'LineWidth', 1.5)
ylabel('integral')
grid on

subplot(3,1,3)
plot(eq29.l, eq29.absGradS2, 'LineWidth', 1.5)
ylabel('|grad S|^2')
xlabel('l')
grid on

figure;
subplot(3,1,1)
plot(eq31.l, eq31.dQdrho, 'LineWidth', 1.5)
ylabel('dQ/drho')
grid on

subplot(3,1,2)
plot(eq31.l, eq31.curvatureTerm, 'LineWidth', 1.5)
ylabel('Eq. (31)')
grid on

subplot(3,1,3)
plot(eq31.l, eq31.driveTerm, 'LineWidth', 1.5)
ylabel('2 \mu_0 p'' Eq.(31)')
xlabel('l')
grid on

mats = assemble_eq32_matrices(geom, eq29, eq31);

sol = solve_eq32_eig(mats, 'NumModes', 8);

figure;
for k = 1:min(4, size(sol.V,2))
    subplot(4,1,k)
    plot(sol.l, sol.V(:,k), 'LineWidth', 1.5);
    title(sprintf('Mode %d, lambda = %.6e', k, real(sol.lambda(k))));
    grid on;
end
xlabel('l');