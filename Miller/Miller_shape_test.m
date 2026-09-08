clc
clear
close all

%% Read GEQDSK and select flux surfaces
%eq = read_geqdsk('./geqdsk_PT0.7');
%eq = read_geqdsk('./geqdsk_circular');
eq = read_geqdsk('./geqdsk_NT0.6');

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
psiN = 0.5;
param = fit_Miller(eq, psiN, 'dpsi', 1.e-3, 'NTheta', 600);
ntheta = 300;

figure;
s = extract_flux_surface(eq, psiN);
[R, Z, u, p] = DshapeParam(param, ntheta);
plot(s.R, s.Z, '-k', 'LineWidth', 2);
hold on;
plot(R, Z, '-r', 'LineWidth', 2);
plot(eq.rbbbs, eq.zbbbs, '-b', 'LineWidth', 1);
xlabel('R [m]');
ylabel('Z [m]');
legend('Fitted surface', 'geqdsk surface')

%% Evaluate poloidal field
[Bp, out] = Bpol_Dshape(param, ntheta);
Bp_EFIT = Bpol_from_EFIT(eq, out.R, out.Z);

figure;
plot(out.u, Bp, '-k', 'LineWidth', 2);
hold on
plot(out.u, Bp_EFIT, '-r', 'LineWidth', 2);
legend('From fitting', 'From geqdsk');

%% Evaluate PEST coordinate theta
integrand = 1./(out.R.^2 .* Bp);
theta_PEST = zeros(1, ntheta);
for i = 2:ntheta
    theta_PEST(i) = theta_PEST(i-1) + 0.5*hypot(out.R(i-1) - out.R(i), ...
        out.Z(i-1) - out.Z(i))*(integrand(i-1) + integrand(i));
end
theta_PEST = 2.*pi.*theta_PEST./theta_PEST(end);