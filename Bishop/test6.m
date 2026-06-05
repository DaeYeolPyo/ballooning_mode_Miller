clc
clear
close all

%% Ballooning stability from Bishop Chapter 3 terms
eq = read_geqdsk('./geqdsk_NT0.6');
%eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
psiN = 0.9;
NTheta = 512;

geom = build_bishop_geometry(eq, psiN, ...
    'NTheta', NTheta, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'UseUniformArclength', true);

%% Single stability solve on a finite ballooning-theta window
nside = 5;
ch3 = compute_bishop_ch3_terms(geom, ...
    'L0Index', 1, ...
    'NPeriodsEachSide', nside, ...
    'UseMu0Pprime', true);

stab = solve_bishop_ballooning_stability(ch3, ...
    'NumModes', 6, ...
    'UseExtended', true, ...
    'Plot', true);

fprintf('Ballooning stability result\n');
fprintf('  psiN = %.4f\n', psiN);
fprintf('  theta window = [%.6g, %.6g]\n', ...
    min(stab.theta), max(stab.theta));
fprintf('  lambda_min = %.12e\n', stab.lambdaMin);
fprintf('  status = %s\n', stab.status);
fprintf('  deltaW lowest = %.12e\n', stab.energy.deltaW(1));
fprintf('  bending lowest = %.12e\n', stab.energy.bending(1));
fprintf('  pressure-curvature lowest = %.12e\n', ...
    stab.energy.pressureCurvature(1));

%% Window-size convergence check
nsideList = 2:7;
lambdaMin = zeros(size(nsideList));
status = strings(size(nsideList));

for i = 1:numel(nsideList)
    ch3i = compute_bishop_ch3_terms(geom, ...
        'L0Index', 1, ...
        'NPeriodsEachSide', nsideList(i), ...
        'UseMu0Pprime', true);

    stabi = solve_bishop_ballooning_stability(ch3i, ...
        'NumModes', 3, ...
        'UseExtended', true, ...
        'Plot', false);

    lambdaMin(i) = stabi.lambdaMin;
    status(i) = string(stabi.status);

    fprintf('  nside=%d: lambda_min=%.12e, status=%s\n', ...
        nsideList(i), lambdaMin(i), status(i));
end

figure;
plot(nsideList, lambdaMin, 'o-', 'LineWidth', 1.4);
yline(0, 'k--');
grid on;
xlabel('NPeriodsEachSide');
ylabel('lambda_min');
title('Ballooning theta window convergence');
