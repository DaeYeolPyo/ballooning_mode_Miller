%% verify_Solovev.m
% End-to-end verification of solve_GS against an analytic Solov'ev solution.

clearvars;
close all;
clc;

%% Exact Solov'ev equilibrium
R0 = 4.0;
B0 = 2.5;
a = 1.0;
kappa = 1.5;
delta = 0.4;
targetPsiAxis = 1.5;

solov = make_Solovev_case(R0, a, kappa, delta, targetPsiAxis);

%% Generating inputs
nBoundary = 1000;
theta = (0:nBoundary-1).'*(2*pi/nBoundary);
boundary = solov.boundary(theta);

nProfile = 101;
profiles = solov.makeProfiles(nProfile);

NR = 128; NZ = 128;
input.grid = struct( ...
    'nBoundary', nBoundary, ...
    'NR', 128, ...
    'NZ', 128, ...
    'quadratureDegree', int8(10), ...
    'nProfile', nProfile);

input.dim.R0 = R0;
input.dim.B0 = B0;
input.dim.aMinor = a;
input.profile = struct( ...
    'mode', 1, ...
    'pprime', profiles.dp_dpsiN, ...
    'FFprime', profiles.FdF_dpsiN, ...
    'pedge', 0.0);

input.boundary = struct( ...
    'R', boundary(:,1), ...
    'Z', boundary(:,2), ...
    'psiBoundary', solov.psiBoundary);

input.Picard = struct( ...
    'omega', 0.5, ...
    'maxIterations', 100, ...
    'updateTolerance', 1.e-9, ...
    'residualTolerance', 1.e-9, ...
    'axisMode', solov.axisMode, ...
    'symmetryTolerance', 1.e-12, ...
    'fluxSpanTolerance', 1.e-12);

input.output = struct( ...
    'verbose', true, ...
    'checkTime', true, ...
    'showPlot', true);

%% GS solve
equilibrium = solve_GS(input);

%% Exact Solov'ev equilibrium
points = equilibrium.mesh.points;
psiExact = solov.psi(points(:,1), points(:,2));
psiError = equilibrium.psi-psiExact;

%% Output report
relativeNodalL2 = norm(psiError)/max(norm(psiExact), eps);
relativeInfinity = norm(psiError, inf)/max(norm(psiExact, inf), eps);
boundaryNodes = equilibrium.mesh.boundaryNodes;
boundaryConditionError = norm( ...
    equilibrium.psi(boundaryNodes)-solov.psiBoundary, inf);
geometryBoundaryMismatch = norm( ...
    psiExact(boundaryNodes)-solov.psiBoundary, inf);
relativeGeometryBoundaryMismatch = ...
    geometryBoundaryMismatch/max(abs(solov.dpsi), eps);
relativeAxisFluxError = abs(equilibrium.solverResult.psiAxis-solov.psiAxis) ...
    /max(abs(solov.psiAxis), eps);
axisPositionError = norm(equilibrium.axisPoint-solov.axisPoint);

%% Safety-factor profile verification
qCheckPsiN = linspace(0,1,21).';
qNumerical = interp1( ...
    equilibrium.profile.psiN, equilibrium.profile.q, ...
    qCheckPsiN, 'pchip');
qExact = exact_Solovev_q( ...
    solov, R0*B0, qCheckPsiN, boundary);
qError = qNumerical-qExact;

relativeQL2 = norm(qError)/max(norm(qExact),eps);
relativeQInfinity = norm(qError,inf)/max(norm(qExact,inf),eps);
relativeQAxisError = abs(qNumerical(1)-qExact(1)) ...
    /max(abs(qExact(1)),eps);
relativeQBoundaryError = abs(qNumerical(end)-qExact(end)) ...
    /max(abs(qExact(end)),eps);

fprintf('\nSolov''ev verification errors\n');
fprintf('  Nodal relative L2       = %.3e\n', relativeNodalL2);
fprintf('  Relative infinity       = %.3e\n', relativeInfinity);
fprintf('  Boundary-condition error= %.3e\n', boundaryConditionError);
fprintf('  Polygon boundary mismatch= %.3e\n', geometryBoundaryMismatch);
fprintf('  Relative axis-flux      = %.3e\n', relativeAxisFluxError);
fprintf('  Axis-position error     = %.3e\n', axisPositionError);
fprintf('  q relative L2           = %.3e\n', relativeQL2);
fprintf('  q relative infinity     = %.3e\n', relativeQInfinity);
fprintf('  Relative q-axis         = %.3e\n', relativeQAxisError);
fprintf('  Relative q-boundary     = %.3e\n', relativeQBoundaryError);

assert(equilibrium.solverResult.converged, ...
    'Picard iteration did not converge.');
assert(equilibrium.solverResult.residualHistory(end) ...
    < input.Picard.residualTolerance, 'Final residual is too large.');
assert(boundaryConditionError < 1.e-13, ...
    'Boundary condition is violated.');
assert(relativeGeometryBoundaryMismatch < 5.e-5, ...
    'Polygon boundary mismatch is too large.');
assert(relativeNodalL2 < 1.e-3, 'Nodal relative L2 error is too large.');
assert(relativeInfinity < 5.e-3, 'Relative infinity error is too large.');
assert(relativeAxisFluxError < 5.e-3, ...
    'Relative axis-flux error is too large.');
assert(axisPositionError < 5.e-2, 'Axis-position error is too large.');
assert(relativeQL2 < 1.e-3, 'q relative L2 error is too large.');
assert(relativeQInfinity < 2.e-3, ...
    'q relative infinity error is too large.');
assert(relativeQAxisError < 2.e-3, ...
    'Relative q-axis error is too large.');
assert(relativeQBoundaryError < 2.e-3, ...
    'Relative q-boundary error is too large.');

plot_Solovev_verification( ...
    equilibrium, solov, boundary, ...
    qCheckPsiN, qNumerical, qExact);

fprintf('\nAll solve_GS Solov''ev verification checks passed.\n');

%% Auxiliary functions
function solov = make_Solovev_case(R0, a, kappa, delta, targetPsiAxis)
    Rinner = R0-a;
    Router = R0+a;
    Rtop = R0-delta*a;
    Ztop = kappa*a;

    A1 = -2/(R0^2*a^2);
    A2 = -2/(kappa^2*a^2);

    h = @(R,Z) [1, R.^2, R.^4-4*R.^2.*Z.^2, R.^2.*log(R)-Z.^2];
    dhdR = @(R,Z) [0, 2*R, 4*R.^3-8*R.*Z.^2, 2*R.*log(R)+R];
    psiPart = @(R,Z) (A1/8).*R.^4+(A2/2).*Z.^2;
    dpsiPartdR = @(R,Z) (A1/2).*R.^3+0.*Z;

    M = [h(Rinner,0); h(Router,0); h(Rtop,Ztop); dhdR(Rtop,Ztop)];
    c = M\-[psiPart(Rinner,0); psiPart(Router,0); ...
        psiPart(Rtop,Ztop); dpsiPartdR(Rtop,Ztop)];

    psiRaw = @(R,Z) (A1/8).*R.^4+(A2/2).*Z.^2+c(1)+c(2).*R.^2 ...
        +c(3).*(R.^4-4*R.^2.*Z.^2)+c(4).*(R.^2.*log(R)-Z.^2);
    dpsiRawdR = @(R,Z) (A1/2).*R.^3+2*c(2).*R ...
        +c(3).*(4*R.^3-8*R.*Z.^2)+c(4).*(2*R.*log(R)+R)+0.*Z;

    Raxis = fzero(@(R) dpsiRawdR(R,0), [Rinner, Router]);
    scale = targetPsiAxis/psiRaw(Raxis,0);

    params = struct( ...
        'A1', scale*A1, 'A2', scale*A2, ...
        'c1', scale*c(1), 'c2', scale*c(2), ...
        'c3', scale*c(3), 'c4', scale*c(4), ...
        'axisBracket', [Rinner, Router], ...
        'rhoMax', 2.5*max(a,kappa*a), ...
        'psiBoundary', 0.0, ...
        'boundaryScanPoints', 500);

    solov = exact_Solovev(params);
end


function plot_Solovev_verification( ...
    equilibrium, solov, boundary, ...
    qCheckPsiN, qNumerical, qExact)
    P3 = equilibrium.mesh;
    P1 = triangulation(P3.P1elements, P3.P1points);

    Rvec = linspace(min(boundary(:,1)), max(boundary(:,1)), 180);
    Zvec = linspace(min(boundary(:,2)), max(boundary(:,2)), 240);
    [Rplot,Zplot] = meshgrid(Rvec,Zvec);

    psiNumerical = eval_P3_sol( ...
        P1, P3, equilibrium.psi, [Rplot(:),Zplot(:)]);
    psiNumerical = reshape(psiNumerical, size(Rplot));

    psiNnumerical = (psiNumerical-equilibrium.solverResult.psiAxis) ...
        /equilibrium.solverResult.dpsi;
    psiNexact = solov.psiN(Rplot,Zplot);
    normalizedError = (psiNumerical-solov.psi(Rplot,Zplot)) ...
        /max(abs(solov.psiAxis), eps);

    outside = ~isfinite(psiNumerical);
    psiNexact(outside) = NaN;
    normalizedError(outside) = NaN;

    Rlimits = [min(Rvec),max(Rvec)];
    Zlimits = [min(Zvec),max(Zvec)];
    fluxLevels = 0:0.1:1;

    figure('Name','Solov''ev field comparison', 'NumberTitle','off', ...
        'Color','w', 'Units','normalized', ...
        'Position',[0.02,0.08,0.96,0.80]);
    layout = tiledlayout(1,3, ...
        'TileSpacing','compact', 'Padding','compact');

    ax1 = nexttile(layout);
    contourf(ax1, Rplot, Zplot, psiNnumerical, fluxLevels, ...
        'LineStyle','none');
    hold(ax1,'on');
    plot(ax1, boundary(:,1), boundary(:,2), 'k-', 'LineWidth',1.5);
    plot(ax1, equilibrium.axisPoint(1), equilibrium.axisPoint(2), ...
        'rx', 'MarkerSize',10, 'LineWidth',2);
    title(ax1,'Numerical \psi_N');
    clim(ax1,[0,1]);
    colorbar(ax1);

    ax2 = nexttile(layout);
    contourf(ax2, Rplot, Zplot, psiNexact, fluxLevels, ...
        'LineStyle','none');
    hold(ax2,'on');
    plot(ax2, boundary(:,1), boundary(:,2), 'k-', 'LineWidth',1.5);
    plot(ax2, solov.axisPoint(1), solov.axisPoint(2), ...
        'rx', 'MarkerSize',10, 'LineWidth',2);
    title(ax2,'Exact \psi_N');
    clim(ax2,[0,1]);
    colorbar(ax2);

    ax3 = nexttile(layout);
    contourf(ax3, Rplot, Zplot, normalizedError, 31, ...
        'LineStyle','none');
    hold(ax3,'on');
    plot(ax3, boundary(:,1), boundary(:,2), 'k-', 'LineWidth',1.5);
    finiteError = abs(normalizedError(isfinite(normalizedError)));
    errorLimit = max(finiteError);
    if errorLimit > 0
        clim(ax3,[-errorLimit,errorLimit]);
    end
    title(ax3,'(\psi_h-\psi_{exact})/|\psi_{axis}|');
    colorbar(ax3);

    for ax = [ax1,ax2,ax3]
        axis(ax,'equal');
        xlim(ax,Rlimits);
        ylim(ax,Zlimits);
        grid(ax,'on');
        box(ax,'on');
        xlabel(ax,'R [m]');
        ylabel(ax,'Z [m]');
        ax.FontSize = 11;
    end

    title(layout,'Solov''ev equilibrium verification');

    result = equilibrium.solverResult;
    figure('Name','Solov''ev Picard convergence', 'NumberTitle','off', ...
        'Color','w', 'Units','normalized', ...
        'Position',[0.20,0.16,0.60,0.58]);
    iteration = 1:result.iterations;
    semilogy(iteration, result.updateHistory, 'o-', 'LineWidth',1.3);
    hold on;
    semilogy(iteration, result.residualHistory, 's-', 'LineWidth',1.3);
    grid on;
    box on;
    xlabel('Picard iteration');
    ylabel('Relative error');
    legend('solution update','nonlinear residual', 'Location','best');
    title('Grad-Shafranov Picard convergence');
    set(gca,'FontSize',11);

    figure('Name','Solov''ev q-profile verification', ...
        'NumberTitle','off', 'Color','w', ...
        'Units','normalized', ...
        'Position',[0.14,0.16,0.72,0.58]);
    qLayout = tiledlayout(1,2, ...
        'TileSpacing','compact', 'Padding','compact');

    axQ = nexttile(qLayout);
    plot(axQ, equilibrium.profile.psiN, equilibrium.profile.q, ...
        'b-', 'LineWidth',1.5);
    hold(axQ,'on');
    plot(axQ, qCheckPsiN, qExact, ...
        'ko', 'MarkerSize',5, 'LineWidth',1.1);
    grid(axQ,'on');
    box(axQ,'on');
    xlim(axQ,[0,1]);
    xlabel(axQ,'\psi_N');
    ylabel(axQ,'q');
    legend(axQ,'P3 GS','exact Solov''ev', 'Location','best');
    title(axQ,'Safety-factor profile');
    axQ.FontSize = 11;

    axQE = nexttile(qLayout);
    relativeQError = abs(qNumerical-qExact) ...
        ./max(abs(qExact),eps);
    semilogy(axQE, qCheckPsiN, relativeQError, ...
        'o-', 'LineWidth',1.3, 'MarkerSize',5);
    grid(axQE,'on');
    box(axQE,'on');
    xlim(axQE,[0,1]);
    xlabel(axQE,'\psi_N');
    ylabel(axQE,'relative error');
    title(axQE,'q-profile error');
    axQE.FontSize = 11;

    title(qLayout,'Solov''ev safety-factor verification');
end


function q = exact_Solovev_q(solov, Fboundary, psiNLevels, boundary)
%EXACT_SOLOVEV_Q Reference q from analytic Solov'ev derivatives.

    psiNLevels = psiNLevels(:);
    nLevel = numel(psiNLevels);

    if nLevel < 3 || any(diff(psiNLevels) <= 0) || ...
            abs(psiNLevels(1)) > 10*eps || ...
            abs(psiNLevels(end)-1) > 10*eps
        error(['[VERIFY_SOLOVEV] q-check levels must be strictly ', ...
               'increasing and include psiN = 0 and 1.']);
    end

    F2 = Fboundary^2 ...
        + 2*solov.FFprime.*(psiNLevels-1);

    if any(F2 <= 0)
        error('[VERIFY_SOLOVEV] Exact F^2 is not positive.');
    end

    F = sign(Fboundary)*sqrt(F2);
    q = nan(nLevel,1);

    Raxis = solov.axisPoint(1);
    p = solov.parameters;

    psiRR = ...
          (3*p.A1/2)*Raxis^2 ...
        + 2*p.c2 ...
        + 12*p.c3*Raxis^2 ...
        + p.c4*(2*log(Raxis)+3);

    psiZZ = ...
          p.A2 ...
        - 8*p.c3*Raxis^2 ...
        - 2*p.c4;

    psiRZ = 0.0;
    exactHessian = [psiRR,psiRZ; psiRZ,psiZZ];

    if det(exactHessian) <= 0
        error('[VERIFY_SOLOVEV] Exact axis Hessian is not elliptic.');
    end

    q(1) = F(1)/(Raxis*sqrt(det(exactHessian)));

    nTheta = size(boundary,1);
    theta = (0:nTheta-1).'*(2*pi/nTheta);
    direction = [cos(theta),sin(theta)];
    boundaryRadius = hypot( ...
        boundary(:,1)-Raxis, ...
        boundary(:,2)-solov.axisPoint(2));

    for j = 2:nLevel
        psiNj = psiNLevels(j);
        targetPsi = solov.psiAxis+solov.dpsi*psiNj;
        surfaceRadius = zeros(nTheta,1);

        if abs(psiNj-1) <= 10*eps
            surfaceRadius = boundaryRadius;
        else
            for k = 1:nTheta
                radialEquation = @(rho) solov.psi( ...
                    Raxis+rho*direction(k,1), ...
                    solov.axisPoint(2)+rho*direction(k,2)) ...
                    - targetPsi;

                surfaceRadius(k) = fzero( ...
                    radialEquation,[0,boundaryRadius(k)]);
            end
        end

        surface = solov.axisPoint ...
            + surfaceRadius.*direction;
        surface = [surface;surface(1,:)];

        segment = diff(surface,1,1);
        dl = hypot(segment(:,1),segment(:,2));
        midpoint = 0.5*(surface(1:end-1,:)+surface(2:end,:));

        gradPsiMagnitude = hypot( ...
            solov.dpsi_dR(midpoint(:,1),midpoint(:,2)), ...
            solov.dpsi_dZ(midpoint(:,1),midpoint(:,2)));

        fluxIntegral = sum( ...
            dl./(midpoint(:,1).*gradPsiMagnitude));

        q(j) = F(j)*fluxIntegral/(2*pi);
    end
end
