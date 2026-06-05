clc
clear
close all

%% Bishop s-alpha diagram example
eq = read_geqdsk('./geqdsk_NT0.6');
%eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
psiN = 0.5;

geom = build_bishop_geometry(eq, psiN, ...
    'NTheta', 256, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'UseUniformArclength', true);

% Direct s scan. Internally Iprime is chosen at each (s,alpha) so that the
% Bishop Eq. (29) secular increment Jperiod matches the requested shear.
sValues = linspace(-2, 2, 5);

sa = scan_bishop_s_alpha_diagram(geom, ...
    'ShearValues', sValues, ...
    'AlphaValues', linspace(0, 3, 13), ...
    'NPeriodsEachSide', 3, ...
    'NumModes', 3, ...
    'RefineCritical', true, ...
    'MaxRefineIter', 12, ...
    'Plot', true, ...
    'Verbose', true);

fprintf('s-alpha scan summary\n');
fprintf('  q0 = %.12e\n', sa.q0);
fprintf('  Cq = %.12e\n', sa.Cq);
fprintf('  shearFactor = %.12e\n', sa.shearFactor);
fprintf('  alphaFactor = %.12e\n', sa.alphaFactor);

for i = 1:numel(sa.s)
    fprintf('  s=% .6e: nRoots=%d', sa.s(i), sa.nRoots(i));
    if sa.nRoots(i) > 0
        fprintf(', first alpha=% .12e', sa.firstAlpha(i));
    end
    fprintf('\n');
end
