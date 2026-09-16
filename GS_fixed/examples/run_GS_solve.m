%% run_GS_equilibrium.m
%
% Driver script for a nonlinear fixed-boundary
% Grad-Shafranov equilibrium.
%
% Current solver conventions:
%
%   psi_boundary = 0
%
%   psi_N = (psi-psi_axis)/(psi_boundary-psi_axis)
%
% Profile inputs:
%
%   pprime  = dp/dpsi_N
%   FFprime = F*dF/dpsi_N
%
% The values used below are illustrative physical profiles.

clearvars;
close all;
clc;

totalTimer = tic;

%% ============================================================
% 1. Physical and geometric inputs
% =============================================================

mu0 = 4*pi*1.e-7;

% Miller-type fixed boundary
R0 = 3.0;        % geometric major radius [m]
a = 0.8;         % horizontal minor radius [m]
kappa = 1.5;     % elongation
delta = 0.25;    % triangularity, |delta| < 1

% Vacuum toroidal field at R0
B0 = 3.0;        % [T]

% Current solver requires this value.
psiBoundary = 0.0;

if R0 <= a
    error('Require R0 > a so that the plasma remains at R > 0.');
end

if abs(delta) >= 1
    error('Miller triangularity must satisfy |delta| < 1.');
end

%% ============================================================
% 2. Mesh parameters
% =============================================================

nBoundary = 200;

NR = 55;
NZ = 65;

quadratureDegree = 10;

%% ============================================================
% 3. Pressure and F profiles
%
% psi_N = 0 : magnetic axis
% psi_N = 1 : plasma boundary
%
% Pressure:
%
%   p(psi_N) = pAxis*(1-psi_N)^alphaP
%
% Toroidal-field function:
%
%   F^2(psi_N)
%       = F_boundary^2
%         - deltaF2*(1-psi_N)^alphaF
%
% Therefore:
%
%   dp/dpsi_N
%       = -alphaP*pAxis*(1-psi_N)^(alphaP-1)
%
%   F*dF/dpsi_N
%       = 0.5*d(F^2)/dpsi_N
% =============================================================

nProfile = 201;

psiNProfile = ...
    linspace(0,1,nProfile).';

% Pressure profile
pAxis = 5.0e4;       % [Pa]
alphaP = 2.0;

pressureProfile = ...
    pAxis*(1-psiNProfile).^alphaP;

pprime = ...
    -alphaP*pAxis ...
    *(1-psiNProfile).^(alphaP-1);

% F profile
Fboundary = R0*B0;   % [T m]

% F_boundary^2-F_axis^2
deltaF2 = 0.10;      % [T^2 m^2]
alphaF = 2.0;

F2Profile = ...
    Fboundary^2 ...
    - deltaF2*(1-psiNProfile).^alphaF;

if any(F2Profile <= 0)
    error('The prescribed F^2 profile is nonpositive.');
end

Fprofile = sqrt(F2Profile);

FFprime = ...
    0.5*deltaF2*alphaF ...
    *(1-psiNProfile).^(alphaF-1);

if any(~isfinite(pprime)) || any(~isfinite(FFprime))
    error('The input profiles contain NaN or Inf.');
end

if numel(pprime) ~= numel(FFprime)
    error('pprime and FFprime must have equal lengths.');
end

fprintf('\n');
fprintf('========================================================\n');
fprintf(' Fixed-boundary Grad-Shafranov equilibrium\n');
fprintf('========================================================\n');

fprintf('R0                  = %.6f m\n',R0);
fprintf('a                   = %.6f m\n',a);
fprintf('kappa               = %.6f\n',kappa);
fprintf('delta               = %.6f\n',delta);
fprintf('B0                  = %.6f T\n',B0);
fprintf('p_axis              = %.6e Pa\n',pAxis);
fprintf('F_axis              = %.6e T m\n',Fprofile(1));
fprintf('F_boundary          = %.6e T m\n',Fprofile(end));

%% ============================================================
% 4. Construct fixed Miller boundary
%
%   R(theta) =
%       R0 + a*cos(theta+asin(delta)*sin(theta))
%
%   Z(theta) =
%       kappa*a*sin(theta)
% =============================================================

theta = ...
    (0:nBoundary-1).'*(2*pi/nBoundary);

triangularityAngle = asin(delta);

Rb = R0+a*cos( ...
    theta+triangularityAngle*sin(theta));

Zb = kappa*a*sin(theta);

plasmaBoundary = [Rb,Zb];

if any(Rb <= 0)
    error('The prescribed LCFS reaches R <= 0.');
end

boundaryConstraints = [ ...
    (1:nBoundary).', ...
    [2:nBoundary,1].'];

%% ============================================================
% 5. Generate interior P1 points
% =============================================================

Rgrid = ...
    linspace(min(Rb),max(Rb),NR);

Zgrid = ...
    linspace(min(Zb),max(Zb),NZ);

[RR,ZZ] = meshgrid(Rgrid,Zgrid);

[insideDomain,onBoundary] = ...
    inpolygon(RR,ZZ,Rb,Zb);

strictlyInside = ...
    insideDomain & ~onBoundary;

interiorPoints = [ ...
    RR(strictlyInside), ...
    ZZ(strictlyInside)];

% Include the geometric center explicitly.
centerPoint = [R0,0.0];

distanceFromCenter = hypot( ...
    interiorPoints(:,1)-centerPoint(1), ...
    interiorPoints(:,2)-centerPoint(2));

coordinateScale = max([ ...
    max(Rb)-min(Rb), ...
    max(Zb)-min(Zb), ...
    1]);

centerTolerance = ...
    1.e-12*coordinateScale;

interiorPoints( ...
    distanceFromCenter < centerTolerance,:) = [];

P1points = [ ...
    plasmaBoundary;
    interiorPoints;
    centerPoint];

%% ============================================================
% 6. Generate constrained P1 mesh
% =============================================================

meshTimer = tic;

meshP1 = generate_P1_mesh( ...
    P1points,boundaryConstraints,false);

triangleCenters = incenter(meshP1);

triangleInside = inpolygon( ...
    triangleCenters(:,1), ...
    triangleCenters(:,2), ...
    Rb,Zb);

if any(~triangleInside)
    error([ ...
        'The constrained triangulation contains triangles ', ...
        'outside the prescribed LCFS.']);
end

P1meshTime = toc(meshTimer);

%% ============================================================
% 7. Upgrade to P3 solution mesh
% =============================================================

P3timer = tic;

P3 = generate_P3_nodes(meshP1);

P3meshTime = toc(P3timer);

Ndof = size(P3.points,1);
Nt = size(P3.elements,1);

fprintf('\n');
fprintf('P1 vertices          = %d\n', ...
    size(P3.P1points,1));

fprintf('P1 triangles         = %d\n',Nt);

fprintf('P3 DOFs              = %d\n',Ndof);

fprintf('P3 boundary DOFs     = %d\n', ...
    numel(P3.boundaryNodes));

fprintf('P1 mesh time         = %.3f s\n',P1meshTime);
fprintf('P3 upgrade time      = %.3f s\n',P3meshTime);

assert(all(P3.points(:,1) > 0), ...
    'The P3 mesh contains R <= 0.');

%% ============================================================
% 8. Precompute FEM quadrature
% =============================================================

quad = precompute_quadrature( ...
    quadratureDegree);

fprintf('Quadrature degree    = %d\n', ...
    quadratureDegree);

fprintf('Quadrature points    = %d\n', ...
    numel(quad.w));

%% ============================================================
% 9. Construct a nonzero initial flux guess
%
% The distance to the discrete LCFS gives a smooth positive
% interior profile.
%
% psiAxisGuess only initializes Picard iteration. It does not
% prescribe the final axis flux.
% =============================================================

psiAxisGuess = 0.5;

distanceToBoundary = ...
    inf(Ndof,1);

for k = 1:nBoundary
    distanceToPoint = hypot( ...
        P3.points(:,1)-Rb(k), ...
        P3.points(:,2)-Zb(k));

    distanceToBoundary = min( ...
        distanceToBoundary,distanceToPoint);
end

maximumDistance = ...
    max(distanceToBoundary);

if maximumDistance <= 0
    error('Failed to construct the initial flux guess.');
end

psi0 = ...
    psiAxisGuess ...
    *distanceToBoundary/maximumDistance;

% Enforce the homogeneous fixed-boundary value exactly.
psi0(P3.boundaryNodes) = ...
    psiBoundary;

%% ============================================================
% 10. Picard iteration options
%
% With the profiles used here, the solution has a positive
% interior maximum and psi_boundary = 0.
% =============================================================

options = struct();

options.omega = 0.5;
options.maxIterations = 100;

options.updateTolerance = 1.e-8;
options.residualTolerance = 1.e-8;

options.axisMode = 'max';

options.symmetryTolerance = 1.e-12;
options.fluxSpanTolerance = 1.e-12;

options.verbose = true;

%% ============================================================
% 11. Solve nonlinear Grad-Shafranov equation
% =============================================================

solveTimer = tic;

[psi,picardResult] = Picard_iteration( ...
    P3,quad, ...
    pprime,FFprime, ...
    psi0,options);

solveTime = toc(solveTimer);

if ~picardResult.converged
    warning([ ...
        'The GS solver did not converge. Inspect the iteration ', ...
        'history before using this equilibrium.']);
end

%% ============================================================
% 12. Construct normalized flux
% =============================================================

psiAxis = ...
    picardResult.psiAxis;

dpsi = ...
    picardResult.dpsi;

psiN = ...
    (psi-psiAxis)/dpsi;

psiNMinimum = min(psiN);
psiNMaximum = max(psiN);

% Clipping is used only for profile visualization.
% It is not used inside the GS solve.
psiNForPlot = ...
    min(max(psiN,0),1);

pressureNodal = interp1( ...
    psiNProfile,pressureProfile, ...
    psiNForPlot,'linear');

F2Nodal = interp1( ...
    psiNProfile,F2Profile, ...
    psiNForPlot,'linear');

FNodal = sqrt(F2Nodal);

BphiNodal = ...
    FNodal./P3.points(:,1);

%% ============================================================
% 13. Locate the nodal magnetic axis
% =============================================================

freeNodes = ...
    picardResult.freeNodes;

switch options.axisMode
    case 'max'
        [psiAxisNodal,localAxisIndex] = ...
            max(psi(freeNodes));

    case 'min'
        [psiAxisNodal,localAxisIndex] = ...
            min(psi(freeNodes));

    otherwise
        error('Nodal axis location requires max or min mode.');
end

axisNode = ...
    freeNodes(localAxisIndex);

axisPoint = ...
    P3.points(axisNode,:);

%% ============================================================
% 14. Element-center magnetic field and toroidal current
%
% Assumed flux convention:
%
%   B_R = -(1/R)*dpsi/dZ
%   B_Z =  (1/R)*dpsi/dR
%   B_phi = F/R
%
%   J_phi =
%       R*dp/dpsi
%       + (F*dF/dpsi)/(mu0*R)
% =============================================================

[Ncenter,dNrefCenter] = ...
    construct_P3_shape_functions(1/3,1/3);

elementCenter = ...
    zeros(Nt,2);

psiCenter = ...
    zeros(Nt,1);

psiNCenter = ...
    zeros(Nt,1);

BRcenter = ...
    zeros(Nt,1);

BZcenter = ...
    zeros(Nt,1);

BphiCenter = ...
    zeros(Nt,1);

JphiCenter = ...
    zeros(Nt,1);

elementArea = ...
    zeros(Nt,1);

totalPlasmaCurrent = 0.0;

L1q = 1-quad.xi-quad.eta;
L2q = quad.xi;
L3q = quad.eta;

for e = 1:Nt
    geometryNodes = ...
        P3.P1elements(e,:);

    Ve = ...
        P3.P1points(geometryNodes,:);

    R1 = Ve(1,1);
    Z1 = Ve(1,2);

    R2 = Ve(2,1);
    Z2 = Ve(2,2);

    R3 = Ve(3,1);
    Z3 = Ve(3,2);

    J = [ ...
        R2-R1, R3-R1;
        Z2-Z1, Z3-Z1];

    detJ = det(J);

    if detJ <= 0
        error( ...
            'Invalid element orientation at element %d.',e);
    end

    elementArea(e) = ...
        0.5*detJ;

    elementCenter(e,:) = ...
        mean(Ve,1);

    Rcenter = ...
        elementCenter(e,1);

    ids = ...
        P3.elements(e,:);

    psiElement = ...
        psi(ids);

    % Reference-to-physical gradient transformation
    dNphysical = ...
        J.'\dNrefCenter;

    gradientPsi = ...
        dNphysical*psiElement;

    dpsi_dR = ...
        gradientPsi(1);

    dpsi_dZ = ...
        gradientPsi(2);

    psiCenter(e) = ...
        Ncenter*psiElement;

    psiNCenter(e) = ...
        (psiCenter(e)-psiAxis)/dpsi;

    BRcenter(e) = ...
        -dpsi_dZ/Rcenter;

    BZcenter(e) = ...
        dpsi_dR/Rcenter;

    psiNProfileCenter = ...
        min(max(psiNCenter(e),0),1);

    Fcenter = interp1( ...
        psiNProfile,Fprofile, ...
        psiNProfileCenter,'linear');

    BphiCenter(e) = ...
        Fcenter/Rcenter;

    dp_dpsiN_center = interp1( ...
        psiNProfile,pprime, ...
        psiNCenter(e),'linear','extrap');

    FdF_dpsiN_center = interp1( ...
        psiNProfile,FFprime, ...
        psiNCenter(e),'linear','extrap');

    dp_dpsi_center = ...
        dp_dpsiN_center/dpsi;

    FdF_dpsi_center = ...
        FdF_dpsiN_center/dpsi;

    JphiCenter(e) = ...
          Rcenter*dp_dpsi_center ...
        + FdF_dpsi_center/(mu0*Rcenter);

    % Higher-order quadrature estimate of total plasma current
    Rq = ...
        L1q*R1+L2q*R2+L3q*R3;

    psiQ = ...
        quad.N*psiElement;

    psiNQ = ...
        (psiQ-psiAxis)/dpsi;

    dp_dpsiN_Q = interp1( ...
        psiNProfile,pprime, ...
        psiNQ,'linear','extrap');

    FdF_dpsiN_Q = interp1( ...
        psiNProfile,FFprime, ...
        psiNQ,'linear','extrap');

    JphiQ = ...
          Rq.*(dp_dpsiN_Q/dpsi) ...
        + (FdF_dpsiN_Q/dpsi) ...
            ./(mu0.*Rq);

    totalPlasmaCurrent = ...
        totalPlasmaCurrent ...
        + sum(quad.w.*detJ.*JphiQ);
end

BpolCenter = hypot( ...
    BRcenter,BZcenter);

BtotalCenter = sqrt( ...
    BRcenter.^2 ...
    + BZcenter.^2 ...
    + BphiCenter.^2);

%% ============================================================
% 15. Final report
% =============================================================

fprintf('\n');
fprintf('========================================================\n');
fprintf(' GS equilibrium result\n');
fprintf('========================================================\n');

fprintf('Converged             = %d\n', ...
    picardResult.converged);

fprintf('Picard iterations     = %d\n', ...
    picardResult.iterations);

fprintf('Final residual        = %.3e\n', ...
    picardResult.residualHistory(end));

fprintf('Final update          = %.3e\n', ...
    picardResult.updateHistory(end));

fprintf('psi_axis              = %+14.7e\n', ...
    psiAxis);

fprintf('dpsi                  = %+14.7e\n', ...
    dpsi);

fprintf('Axis node             = %d\n', ...
    axisNode);

fprintf('Axis location         = (%.8f, %.8f) m\n', ...
    axisPoint(1),axisPoint(2));

fprintf('Nodal psiN range      = [%+.6e, %+.6e]\n', ...
    psiNMinimum,psiNMaximum);

fprintf('Estimated plasma I    = %.6e A\n', ...
    totalPlasmaCurrent);

fprintf('Linear solve time     = %.3f s\n', ...
    solveTime);

fprintf('Total driver time     = %.3f s\n', ...
    toc(totalTimer));

fprintf('========================================================\n');

if psiNMinimum < -1.e-3 || psiNMaximum > 1+1.e-3
    warning([ ...
        'The converged nodal psiN range extends outside [0,1]. ', ...
        'Inspect possible overshoot or an incorrect axis mode.']);
end

%% ============================================================
% 16. Package equilibrium result
% =============================================================

equilibrium = struct();

equilibrium.geometry = struct( ...
    'R0',R0, ...
    'a',a, ...
    'kappa',kappa, ...
    'delta',delta, ...
    'boundary',plasmaBoundary);

equilibrium.mesh = P3;

equilibrium.psi = psi;
equilibrium.psiN = psiN;

equilibrium.psiAxis = psiAxis;
equilibrium.psiBoundary = psiBoundary;
equilibrium.dpsi = dpsi;

equilibrium.axisNode = axisNode;
equilibrium.axisPoint = axisPoint;

equilibrium.profile = struct();

equilibrium.profile.psiN = psiNProfile;

equilibrium.profile.pressure = ...
    pressureProfile;

equilibrium.profile.dp_dpsiN = ...
    pprime;

equilibrium.profile.F = ...
    Fprofile;

equilibrium.profile.FdF_dpsiN = ...
    FFprime;

equilibrium.pressureNodal = ...
    pressureNodal;

equilibrium.FNodal = ...
    FNodal;

equilibrium.BphiNodal = ...
    BphiNodal;

equilibrium.elementCenter = ...
    elementCenter;

equilibrium.BRcenter = ...
    BRcenter;

equilibrium.BZcenter = ...
    BZcenter;

equilibrium.BpolCenter = ...
    BpolCenter;

equilibrium.BphiCenter = ...
    BphiCenter;

equilibrium.BtotalCenter = ...
    BtotalCenter;

equilibrium.JphiCenter = ...
    JphiCenter;

equilibrium.totalPlasmaCurrent = ...
    totalPlasmaCurrent;

equilibrium.solverResult = ...
    picardResult;

%% ============================================================
% 17. Construct P3 visualization connectivity
% =============================================================

localSubtriangles = [ ...
     1,  4,  9;
     4, 10,  9;
     4,  5, 10;
     5,  6, 10;
     5,  2,  6;
     9, 10,  8;
    10,  7,  8;
    10,  6,  7;
     8,  7,  3];

Tplot = zeros(9*Nt,3);

for e = 1:Nt
    rows = ...
        (e-1)*9+(1:9);

    ids = ...
        P3.elements(e,:);

    Tplot(rows,:) = ...
        ids(localSubtriangles);
end

Rnodes = P3.points(:,1);
Znodes = P3.points(:,2);

%% ============================================================
% 18. Plot mesh
% =============================================================

figure('Color','w');

triplot( ...
    P3.P1elements, ...
    P3.P1points(:,1), ...
    P3.P1points(:,2), ...
    'Color',[0.75,0.75,0.75]);

hold on;

plot( ...
    Rb,Zb, ...
    'k-','LineWidth',2);

plot( ...
    axisPoint(1),axisPoint(2), ...
    'ro','MarkerFaceColor','r');

axis equal tight;
grid on;

xlabel('R [m]');
ylabel('Z [m]');

legend( ...
    'P1 mesh', ...
    'prescribed LCFS', ...
    'nodal magnetic axis', ...
    'Location','best');

title('Fixed-boundary computational mesh');

%% ============================================================
% 19. Plot equilibrium fields
% =============================================================

figure('Color','w');

subplot(1,3,1);

trisurf( ...
    Tplot, ...
    Rnodes,Znodes,psi, ...
    'EdgeColor','none');

view(2);
axis equal tight;
colorbar;

xlabel('R [m]');
ylabel('Z [m]');
title('\psi');

subplot(1,3,2);

trisurf( ...
    Tplot, ...
    Rnodes,Znodes,psiN, ...
    'EdgeColor','none');

view(2);
axis equal tight;
colorbar;

xlabel('R [m]');
ylabel('Z [m]');
title('\psi_N');

subplot(1,3,3);

trisurf( ...
    Tplot, ...
    Rnodes,Znodes,pressureNodal, ...
    'EdgeColor','none');

view(2);
axis equal tight;
colorbar;

xlabel('R [m]');
ylabel('Z [m]');
title('Pressure [Pa]');

sgtitle('Grad-Shafranov equilibrium');

%% ============================================================
% 20. Plot Picard convergence
% =============================================================

figure('Color','w');

iterationNumber = ...
    1:picardResult.iterations;

semilogy( ...
    iterationNumber, ...
    picardResult.updateHistory, ...
    'o-','LineWidth',1.5);

hold on;

semilogy( ...
    iterationNumber, ...
    picardResult.residualHistory, ...
    's-','LineWidth',1.5);

grid on;

xlabel('Picard iteration');
ylabel('Relative error');

legend( ...
    'solution update', ...
    'nonlinear residual', ...
    'Location','best');

title('Grad-Shafranov Picard convergence');

%% 21. Contour plot
eq = equilibrium;
P3 = eq.mesh;

% P3 연결 관계와 동일한 P1 geometry triangulation
P1 = triangulation( ...
    P3.P1elements, ...
    P3.P1points);

% Contour용 Cartesian grid
Rvec = linspace( ...
    min(P3.P1points(:,1)), ...
    max(P3.P1points(:,1)), 400);

Zvec = linspace( ...
    min(P3.P1points(:,2)), ...
    max(P3.P1points(:,2)), 500);

[Rplot,Zplot] = meshgrid(Rvec,Zvec);

% 각 격자점에서 P3 FEM solution 평가
psiPlot = eval_P3_sol( ...
    P1, P3, eq.psi, ...
    [Rplot(:),Zplot(:)]);

psiPlot = reshape(psiPlot,size(Rplot));

% Normalized flux
psiNPlot = ...
    (psiPlot-eq.psiAxis)/eq.dpsi;

% Flux-surface contours
figure('Color','w');

levels = 0.1:0.1:0.9;

[C,h] = contour( ...
    Rplot,Zplot,psiNPlot,levels, ...
    'LineWidth',1.3);

clabel(C,h,'FontSize',9);

hold on;

% LCFS: psi_N = 1
boundary = eq.geometry.boundary;

plot( ...
    boundary(:,1),boundary(:,2), ...
    'k-','LineWidth',2);

% Magnetic axis
plot( ...
    eq.axisPoint(1),eq.axisPoint(2), ...
    'ro','MarkerFaceColor','r');

axis equal tight;
grid on;
box on;

xlabel('R [m]');
ylabel('Z [m]');
title('Normalized poloidal-flux contours, \psi_N');

colormap(turbo);
colorbar;
clim([0,1]);