clc
clear
close all

%% Bishop Chapter 3 ballooning equation terms
%eq = read_geqdsk('./curavture_analysis/geqdsk/scan_B2.5_C15_G09_H06.geqdsk');
eq = read_geqdsk('./geqdsk_PT0.6');
psiN = 0.5;
NTheta = 512;

geom = build_bishop_geometry(eq, psiN, ...
    'NTheta', NTheta, ...
    'Direction', 'clockwise', ...
    'NormalDirection', 'bishop', ...
    'UseUniformArclength', true);

ch3 = compute_bishop_ch3_terms(geom, ...
    'L0Index', 1, ...
    'NPeriodsEachSide', 3, ...
    'UseMu0Pprime', true, ...
    'Plot', true);

fprintf('Chapter 3 term check\n');
fprintf('  n=%d, L=%.8g, theta=[%.6g %.6g]\n', ...
    numel(ch3.l), ch3.L, min(ch3.theta), max(ch3.theta));
fprintf('  Jperiod=%.8e\n', ch3.eq29.Jperiod);
fprintf('  gradS2=[%.8e %.8e]\n', ...
    min(ch3.eq29.gradS2), max(ch3.eq29.gradS2));
fprintf('  Eq31 drive=[%.8e %.8e]\n', ...
    min(ch3.eq31.drive), max(ch3.eq31.drive));

if ~isempty(ch3.extended)
    ex = ch3.extended;
    fprintf('  extended theta=[%.6g %.6g], n=%d\n', ...
        min(ex.theta), max(ex.theta), numel(ex.theta));
    fprintf('  extended J=[%.8e %.8e]\n', ...
        min(ex.eq29.J), max(ex.eq29.J));
end

%% Optional: evaluate equation terms for a trial localized F(theta)
if ~isempty(ch3.extended)
    theta = ch3.extended.theta;
    F = exp(-0.25 .* theta.^2);

    ch3F = compute_bishop_ch3_terms(geom, ...
        'L0Index', 1, ...
        'NPeriodsEachSide', 3, ...
        'UseMu0Pprime', true, ...
        'ModeF', F);

    mt = ch3F.modeTerms;

    figure;
    plot(theta, mt.fieldLineBending, 'LineWidth', 1.2); hold on;
    plot(theta, mt.pressureCurvature, 'LineWidth', 1.2);
    plot(theta, mt.residual, 'k--', 'LineWidth', 1.2);
    grid on;
    xlabel('theta');
    ylabel('term value');
    title('Ballooning equation terms for trial F');
    legend('field-line bending', 'pressure-curvature', 'residual');
end
