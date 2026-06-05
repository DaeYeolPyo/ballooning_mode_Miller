clc
clear
close all

%% Critical p' scan for Bishop ballooning stability
eq = read_geqdsk('./geqdsk_NT0.6');
%eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
psiN = 0.5;

geom = build_bishop_geometry(eq, psiN, ...
    'NTheta', 384, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'UseUniformArclength', true);

% The default scanned quantity is pprimeEquation = mu0*dp/dpsi_GEQDSK.
% These GEQDSK files have pprimeEquation < 0, so positive scale follows
% the usual pressure-decreasing direction.
scaleValues = linspace(0, 3, 13);

scan = scan_bishop_critical_pprime(geom, ...
    'ScaleValues', scaleValues, ...
    'NPeriodsEachSide', 5, ...
    'NumModes', 3, ...
    'RefineCritical', true, ...
    'Plot', true, ...
    'Verbose', true);

fprintf('Critical pprime scan summary\n');
fprintf('  base mu0*pprime = %.12e\n', scan.basePprimeEquation);

if isempty(scan.critical)
    fprintf('  No marginal crossing found in scanned range.\n');
    fprintf('  lambda_min range = [%.12e, %.12e]\n', ...
        min(scan.lambdaMin, [], 'omitnan'), max(scan.lambdaMin, [], 'omitnan'));
else
    fprintf('  critical scale = %.12e\n', scan.critical.parameter);
    fprintf('  critical mu0*pprime = %.12e\n', scan.critical.pprimeEquation);
    fprintf('  lambda at critical = %.12e\n', scan.critical.lambdaMin);

    if isfield(geom, 'pprime0') && abs(geom.mu0_pprime0) > 0
        pprimeCrit = geom.pprime0 .* scan.critical.parameter;
        fprintf('  critical raw pprime = %.12e\n', pprimeCrit);
    end
end
