clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, '..', 'Gaur'));

%geqdskFile = fullfile(thisDir, '..', 'Bishop', 'curavture_analysis', ...
%    'geqdsk', 'scan_B2.5_C16_G08_H13.geqdsk');
geqdskFile = "C:\Users\DaeYeolPyo\Desktop\CShape\c-shape-samples\c-shape-samples\kappa5\scan_B2.5_C10_G13_H12_app2.033_gpp0.935_aff0.938_gff1.800.geqdsk";
psin0 = 0.90;

% Coarse first-pass scan. Increase these after the marginal curve location
% is identified.
ns = 30;
nalpha = 30;
ntheta0 = 15;
sGrid = linspace(0.0, 7.0, ns);
alphaGrid = linspace(0.0, 20.0, nalpha);
theta0Grid = linspace(0.0, pi, ntheta0);

eq = read_geqdsk(geqdskFile);
base = cshape_local_equilibrium(eq, psin0, ...
    'NTheta', 161, ...
    'NGeom', 801, ...
    'UseBpFit', true);

fprintf('C-shape target surface diagnostic\n');
fprintf('  file       : %s\n', geqdskFile);
fprintf('  psin       : %.6f\n', psin0);
fprintf('  q profile  : %.8f\n', base.q_profile);
fprintf('  q check    : %.8f\n', base.q_check);
fprintf('  alphaFactor: %.12e\n', base.alphaFactor);
fprintf('  [A B C G H]: [% .8f % .8f % .8f % .8f % .8f]\n', base.coeffs(:,1));
fprintf('dr[A B C G H]: [% .8f % .8f % .8f % .8f % .8f]\n', base.coeffs(:,2));

figure('Color', 'w', 'Name', 'C-shape target surface and Bp');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
contour(eq.rgrid, eq.zgrid, eq.psirz.', 40);
hold on;
surf = extract_flux_surface(eq, psin0);
plot(surf.R, surf.Z, 'k-', 'LineWidth', 1.5);
plot([base.R; base.R(1)], [base.Z; base.Z(1)], 'r--', 'LineWidth', 1.5);
axis equal;
grid on;
xlabel('R');
ylabel('Z');
legend('EFIT contours', 'target contour', 'C-shape fit on PEST \theta', ...
    'Location', 'best');
title('Flux-surface fit');

nexttile;
thetaClosed = [base.theta; 2*pi];
plot(thetaClosed, [base.Bp; base.Bp(1)], 'LineWidth', 1.5);
grid on;
xlabel('\theta_{PEST}');
ylabel('B_p');
title('Fitted local poloidal field');

out = scan_cshape_salpha_diagram(eq, psin0, ...
    'SGrid', sGrid, ...
    'AlphaGrid', alphaGrid, ...
    'Theta0Grid', theta0Grid, ...
    'NThetaCoeff', 161, ...
    'NGeom', 801, ...
    'ThetaB', 5*pi, ...
    'NBalloon', 201, ...
    'NEigs', 6, ...
    'UseBpFit', true, ...
    'ShearDefinition', 'flux', ...
    'Verbose', true);

plot_cshape_salpha_diagram(out);

fprintf('C-shape s-alpha scan done.\n');
fprintf('lambda_max range: %.6g to %.6g\n', ...
    min(out.lambda_max(:), [], 'omitnan'), ...
    max(out.lambda_max(:), [], 'omitnan'));
