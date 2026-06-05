clc
clear
close all

%% Stable/unstable classification on a Bishop s-alpha grid
eq = read_geqdsk('./geqdsk_PT0.6');
%eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
psiN = 0.9;

geom = build_bishop_geometry(eq, psiN, ...
    'NTheta', 256, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'UseUniformArclength', true);

sValues = linspace(0.0, 2.0, 50);
alphaValues = linspace(0, 4.0, 50);

grid = classify_bishop_s_alpha_grid(geom, ...
    'SValues', sValues, ...
    'AlphaValues', alphaValues, ...
    'NPeriodsEachSide', 3, ...
    'NumModes', 3, ...
    'Plot', true, ...
    'Verbose', true);

fprintf('s-alpha classification summary\n');
fprintf('  s definition: s = shearFactor*(dq/dpsi)/q0\n');
fprintf('  shearFactor = %.12e\n', grid.constants.shearFactor);
fprintf('  alpha definition: alpha = alphaFactor*(mu0*dp/dpsi_GEQDSK)\n');
fprintf('  alphaFactor = %.12e\n', grid.constants.alphaFactor);

for i = 1:numel(grid.s)
    for j = 1:numel(grid.alpha)
        fprintf('  s=% .4f, alpha=% .4f -> %s (lambda=%.6e)\n', ...
            grid.s(i), grid.alpha(j), grid.status(i,j), grid.lambdaMin(i,j));
    end
end
