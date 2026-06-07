clc
clear
close all

%% Load ballooning solver
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, '..', 'Gaur'));

%% Read GEQDSK and fit the flux surface
eq = read_geqdsk('./geqdsk_PT0.7');

psiN = 0.77;
p = fit_Miller(eq, psiN, 'dpsi', 1.e-3, 'NTheta', 600);
ntheta = 300;

% Coarse default scan. Increase these grids after the marginal curve appears.
ns = 20; nalpha = 20; ntheta0 = 15;
sGrid = linspace(0.0, 7.0, ns);
alphaGrid = linspace(0.0, 10.0, nalpha);
theta0Grid = linspace(0.0, pi, ntheta0);

out = scan_miller_salpha_diagram(p, ...
    'SGrid', sGrid, ...
    'AlphaGrid', alphaGrid, ...
    'Theta0Grid', theta0Grid, ...
    'NThetaCoeff', 161, ...
    'NGeom', 601, ...
    'ThetaB', 5*pi, ...
    'NBalloon', 201, ...
    'NEigs', 6, ...
    'Verbose', true);

plot_salpha_diagram(out, 'MarginalMode', 'curves');

fprintf('Miller D-shape s-alpha scan done.\n');
fprintf('lambda_max range: %.6g to %.6g\n', ...
    min(out.lambda_max(:), [], 'omitnan'), ...
    max(out.lambda_max(:), [], 'omitnan'));
