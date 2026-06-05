%RUN_MILLER_DSHAPE_SALPHA Validate Bishop ballooning scan on Miller D-shape data.
%
% This script imports the D-shape example from ../Miller/test_Dshape_Miller.m,
% converts its shape and Bp to Bishop's clockwise local coordinate, and
% classifies stable/unstable points on an (s, alpha) grid.

clc;
clear;
close all;

%%
millerPath = fullfile(fileparts(pwd), 'Miller');
nTheta = 128;
sourceNTheta = 1001;
numBallooningPhases = 8;
nPeriodsEachSide = 4;

geom = build_bishop_geometry_from_miller_dshape( ...
    'MillerPath', millerPath, ...
    'SourceNTheta', sourceNTheta, ...
    'NTheta', nTheta, ...
    'JShearFactor', -2*pi);

figure;

tiledlayout(4, 2, "TileSpacing", "tight", ...
    "Padding", "tight");
nexttile(1, [4 1]);
plot(geom.R, geom.Z, 'k-', 'LineWidth', 1.5); hold on;
xlabel('R [m]');
ylabel('Z [m]');

plot(geom.R(1), geom.Z(1), 'mo')

quiver(geom.R, geom.Z, geom.tR, geom.tZ, 0.15, ...
    'b', 'LineWidth', 1.2);
quiver(geom.R, geom.Z, geom.nR, geom.nZ, 0.15, ...
    'r', 'LineWidth', 1.2);

legend('Flux surface', 'First point', 'Tangent', 'Bishop normal');
grid on;

nexttile(2);
plot(geom.l, geom.Bp);
legend('B_p(l)');
grid on;

nexttile(4);
plot(geom.l, geom.R);
hold on;
plot(geom.l, geom.Z);
legend('R(l)', 'Z(l)');
grid on;

nexttile(6);
plot(geom.l, geom.u);
hold on;
plot(geom.l, geom.h0);
plot(geom.l, geom.Rc);
legend('u(l)', 'h_0(l)', 'R_c(l)');
grid on;

p = geom.miller.params;

nexttile(8);
axis off;
txt = sprintf(['A = %f \\kappa = %f\n', ...
    '\\delta = %f\n', ...
    's = %f \\alpha = %f\n'], ...
    p.A, p.kappa, p.delta, p.s_hat, p.alpha);
text(0.05, 0.9, txt, ...
    'Units', 'normalized', ...
    'FontSize', 12, ...
    'VerticalAlignment', 'top', ...
    'Interpreter', 'tex');


%%

% The Miller paper example has s_hat=2.47, so the default range extends
% slightly above the 0 <= s <= 2 range in the reference-style plots.
sValues = linspace(0, 8.0, 41);
alphaValues = linspace(0, 10.0, 41);

grid = classify_bishop_s_alpha_grid(geom, ...
    'SValues', sValues, ...
    'AlphaValues', alphaValues, ...
    'AlphaFactor', geom.alphaFactor0, ...
    'Q0', geom.q0, ...
    'ShearMapping', 'jperiod', ...
    'JShearFactor', geom.JShearFactor0, ...
    'ScanBallooningPhase', true, ...
    'NumBallooningPhases', numBallooningPhases, ...
    'NPeriodsEachSide', nPeriodsEachSide, ...
    'NumModes', 1, ...
    'Plot', false, ...
    'Verbose', false);

resultPath = fullfile(pwd, 'miller_dshape_salpha_result.mat');
save(resultPath, 'geom', 'grid', 'sValues', 'alphaValues');

fig = figure('Color', 'w', 'Name', 'Miller D-shape Bishop s-alpha validation');
tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(geom.R, geom.Z, 'k-', 'LineWidth', 1.8);
hold on;
plot(geom.R(1), geom.Z(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
axis equal;
grid on;
xlabel('R/r');
ylabel('Z/r');
title(sprintf('\\kappa=%.2f, \\delta=%.3f', p.kappa, p.delta));

nexttile;
theta = 2*pi .* geom.l ./ geom.L;
plot(theta, geom.Bp, 'LineWidth', 1.6);
grid on;
xlabel('\theta_l');
ylabel('B_p/B_0');
title(sprintf('q=%.2f, q_{check}=%.2f', p.q, geom.miller.q_check));

nexttile;
imagesc(grid.alpha, grid.s, double(grid.unstable));
set(gca, 'YDir', 'normal');
hold on;
colormap(gca, [0.22 0.43 0.76; 0.82 0.25 0.22]);
caxis([0 1]);
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'stable', 'unstable'};

if any(isfinite(grid.lambdaMin(:)))
    contour(grid.alpha, grid.s, grid.lambdaMin, [0 0], ...
        'k-', 'LineWidth', 1.4);
end

if p.alpha >= min(alphaValues) && p.alpha <= max(alphaValues) && ...
        p.s_hat >= min(sValues) && p.s_hat <= max(sValues)
    plot(p.alpha, p.s_hat, 'kp', ...
        'MarkerFaceColor', 'y', 'MarkerSize', 12, 'LineWidth', 1.2);
end

grid on;
xlabel('\alpha');
ylabel('s');
title('Bishop ballooning stability (min over l_0)');
axis([min(alphaValues) max(alphaValues) min(sValues) max(sValues)]);

pngPath = fullfile(pwd, 'miller_dshape_salpha.png');
try
    exportgraphics(fig, pngPath, 'Resolution', 220);
catch
    saveas(fig, pngPath);
end

fprintf('\nMiller D-shape Bishop validation complete.\n');
fprintf('  A=%.4f, kappa=%.4f, delta=%.4f\n', p.A, p.kappa, p.delta);
fprintf('  q target=%.6f, q check=%.6f\n', p.q, geom.miller.q_check);
fprintf('  s_hat point=%.6f, alpha point=%.6f\n', p.s_hat, p.alpha);
fprintf('  alphaFactor=%.12e, mu0*pprime baseline=%.12e\n', ...
    grid.constants.alphaFactor, geom.mu0_pprime0);
fprintf('  dpsi/dr=%.12e\n', geom.miller.dpdr_flux);
fprintf('  JShearFactor=%.12e\n', grid.constants.JShearFactor);
fprintf('  ballooning phase scan: %d theta_k values\n', ...
    numel(grid.ballooningPhaseIndices));
fprintf('  l0Index range selected: [%g, %g]\n', ...
    min(grid.phaseIndexMin(:), [], 'omitnan'), ...
    max(grid.phaseIndexMin(:), [], 'omitnan'));
fprintf('  result mat: %s\n', resultPath);
fprintf('  result png: %s\n', pngPath);
