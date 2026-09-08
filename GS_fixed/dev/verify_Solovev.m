%% verify_Solovev.m
%
% Verification of the P3 Grad-Shafranov solver using the
% general analytic Solov'ev equilibrium
%
%   psi(R,Z) =
%       A1/8 * R^4
%       + A2/2 * Z^2
%       + c1
%       + c2 * R^2
%       + c3 * (R^4 - 4*R^2*Z^2)
%       + c4 * (R^2*log(R) - Z^2)
%
% This is a full script, not a function.

clear;
close all;
clc;

%% ============================================================
% 1. Desired plasma-shape parameters
% =============================================================

Rgeo = 4.0;
a = 1.0;

kappa = 1.5;
delta = -0.6;

targetPsiAxis = 1.e-2;
psiBoundary = 0.0;

mu0 = 4*pi*1.e-7;

% Desired control points
Rinner = Rgeo-a;
Router = Rgeo+a;

Rtop = Rgeo-delta*a;
Ztop = kappa*a;

%% ============================================================
% 2. Select source coefficients A1 and A2
%
% Only their ratio affects the normalized shape.
% The overall magnitude is adjusted below so that
% psi_axis = targetPsiAxis.
% =============================================================

A1base = ...
    -2/(Rgeo^2*a^2);

A2base = ...
    -2/(kappa^2*a^2);

%% ============================================================
% 3. Determine c1, c2, c3, c4 from shape constraints
%
% Four conditions:
%
%   psi(Rinner,0) = 0
%   psi(Router,0) = 0
%   psi(Rtop,Ztop) = 0
%   dpsi/dR(Rtop,Ztop) = 0
%
% The final condition makes the upper point a horizontal
% tangent point of the LCFS.
% =============================================================

homogeneousBasis = @(R,Z) [ ...
    1, ...
    R.^2, ...
    R.^4-4*R.^2.*Z.^2, ...
    R.^2.*log(R)-Z.^2];

homogeneousBasis_dR = @(R,Z) [ ...
    0, ...
    2*R, ...
    4*R.^3-8*R.*Z.^2, ...
    2*R.*log(R)+R];

particularPsi = @(R,Z) ...
      (A1base/8).*R.^4 ...
    + (A2base/2).*Z.^2;

particularPsi_dR = @(R,Z) ...
    (A1base/2).*R.^3+0.*Z;

shapeMatrix = [ ...
    homogeneousBasis(Rinner,0);
    homogeneousBasis(Router,0);
    homogeneousBasis(Rtop,Ztop);
    homogeneousBasis_dR(Rtop,Ztop)];

shapeRhs = -[ ...
    particularPsi(Rinner,0);
    particularPsi(Router,0);
    particularPsi(Rtop,Ztop);
    particularPsi_dR(Rtop,Ztop)];

conditionNumber = cond(shapeMatrix);

if conditionNumber > 1.e12
    error(['The Solov''ev shape-fitting matrix is nearly ', ...
           'singular. cond(M) = %.3e'],conditionNumber);
end

cBase = shapeMatrix\shapeRhs;

c1base = cBase(1);
c2base = cBase(2);
c3base = cBase(3);
c4base = cBase(4);

shapeConstraintResidual = ...
    norm(shapeMatrix*cBase-shapeRhs,inf);

%% ============================================================
% 4. Find the magnetic axis of the unscaled solution
% =============================================================

psiRaw = @(R,Z) ...
      (A1base/8).*R.^4 ...
    + (A2base/2).*Z.^2 ...
    + c1base ...
    + c2base.*R.^2 ...
    + c3base.*(R.^4-4.*R.^2.*Z.^2) ...
    + c4base.*(R.^2.*log(R)-Z.^2);

dpsiRaw_dR = @(R,Z) ...
      (A1base/2).*R.^3 ...
    + 2*c2base.*R ...
    + c3base.*(4.*R.^3-8.*R.*Z.^2) ...
    + c4base.*(2.*R.*log(R)+R) ...
    + 0.*Z;

axisBracket = [Rinner,Router];

RaxisRaw = fzero( ...
    @(R) dpsiRaw_dR(R,0.0), ...
    axisBracket);

psiAxisRaw = psiRaw(RaxisRaw,0.0);

if psiAxisRaw <= psiBoundary
    error(['The fitted equilibrium does not contain a ', ...
           'positive-flux magnetic axis.']);
end

%% ============================================================
% 5. Scale all coefficients to prescribe psi_axis
%
% Multiplying every coefficient by the same number does not
% change the LCFS geometry.
% =============================================================

coefficientScale = ...
    targetPsiAxis/psiAxisRaw;

A1 = coefficientScale*A1base;
A2 = coefficientScale*A2base;

c1 = coefficientScale*c1base;
c2 = coefficientScale*c2base;
c3 = coefficientScale*c3base;
c4 = coefficientScale*c4base;

%% ============================================================
% 6. Construct analytic Solov'ev model
% =============================================================

params = struct();

params.A1 = A1;
params.A2 = A2;

params.c1 = c1;
params.c2 = c2;
params.c3 = c3;
params.c4 = c4;

params.mu0 = mu0;
params.psiBoundary = psiBoundary;

params.axisBracket = axisBracket;

% Must be larger than the maximum expected distance from
% the magnetic axis to the LCFS.
params.rhoMax = ...
    2.5*max(a,kappa*a);

params.boundaryScanPoints = 500;

solov = exact_Solovev(params);

fprintf('\n');
fprintf('========================================================\n');
fprintf(' General analytic Solov''ev equilibrium\n');
fprintf('========================================================\n');

fprintf('A1                     = %+14.7e\n',solov.A1);
fprintf('A2                     = %+14.7e\n',solov.A2);

fprintf('c1                     = %+14.7e\n',solov.c1);
fprintf('c2                     = %+14.7e\n',solov.c2);
fprintf('c3                     = %+14.7e\n',solov.c3);
fprintf('c4                     = %+14.7e\n',solov.c4);

fprintf('Shape matrix condition = %.3e\n',conditionNumber);
fprintf('Shape fitting residual = %.3e\n', ...
    shapeConstraintResidual);

fprintf('Exact axis             = (%.8f, %.8f)\n', ...
    solov.axisPoint(1),solov.axisPoint(2));

fprintf('psi_axis               = %+14.7e\n', ...
    solov.psiAxis);

fprintf('psi_boundary           = %+14.7e\n', ...
    solov.psiBoundary);

fprintf('dpsi                   = %+14.7e\n', ...
    solov.dpsi);

fprintf('Axis mode              = %s\n', ...
    solov.axisMode);

fprintf('Physical dp/dpsi       = %+14.7e\n', ...
    solov.pprimePhysical);

fprintf('Physical F*dF/dpsi     = %+14.7e\n', ...
    solov.FFprimePhysical);

fprintf('dp/dpsiN               = %+14.7e\n', ...
    solov.pprime);

fprintf('F*dF/dpsiN             = %+14.7e\n', ...
    solov.FFprime);

%% ============================================================
% 7. Check the prescribed shaping conditions
% =============================================================

innerFluxError = ...
    abs(solov.psi(Rinner,0)-psiBoundary);

outerFluxError = ...
    abs(solov.psi(Router,0)-psiBoundary);

topFluxError = ...
    abs(solov.psi(Rtop,Ztop)-psiBoundary);

topTangencyError = ...
    abs(solov.dpsi_dR(Rtop,Ztop));

fprintf('\n');
fprintf('Inner-point flux error = %.3e\n',innerFluxError);
fprintf('Outer-point flux error = %.3e\n',outerFluxError);
fprintf('Top-point flux error   = %.3e\n',topFluxError);
fprintf('Top tangency error     = %.3e\n',topTangencyError);

assert(innerFluxError < 1.e-12);
assert(outerFluxError < 1.e-12);
assert(topFluxError < 1.e-12);
assert(topTangencyError < 1.e-12);

%% ============================================================
% 8. Generate the analytic LCFS
% =============================================================

nBoundary = 180;

theta = ...
    (0:nBoundary-1).'*(2*pi/nBoundary);

boundary = solov.boundary(theta);

Rb = boundary(:,1);
Zb = boundary(:,2);

analyticBoundaryError = ...
    max(abs(solov.psi(Rb,Zb)-psiBoundary));

fprintf('Analytic LCFS error    = %.3e\n', ...
    analyticBoundaryError);

assert(analyticBoundaryError < 1.e-11, ...
    'Analytic boundary does not satisfy psi=psiBoundary.');

%% ============================================================
% 9. Generate interior P1 points
% =============================================================

NR = 55;
NZ = 55;

Rgrid = linspace(min(Rb),max(Rb),NR);
Zgrid = linspace(min(Zb),max(Zb),NZ);

[RR,ZZ] = meshgrid(Rgrid,Zgrid);

[inDomain,onBoundary] = ...
    inpolygon(RR,ZZ,Rb,Zb);

strictlyInside = ...
    inDomain & ~onBoundary;

interiorPoints = [ ...
    RR(strictlyInside), ...
    ZZ(strictlyInside)];

% Explicitly include the exact magnetic axis.
axisPoint = solov.axisPoint;

distanceFromAxis = hypot( ...
    interiorPoints(:,1)-axisPoint(1), ...
    interiorPoints(:,2)-axisPoint(2));

coordinateScale = max([ ...
    max(Rb)-min(Rb), ...
    max(Zb)-min(Zb), ...
    1]);

axisTolerance = ...
    1.e-12*coordinateScale;

interiorPoints( ...
    distanceFromAxis < axisTolerance,:) = [];

%% ============================================================
% 10. Construct the constrained P1 triangulation
% =============================================================

P1points = [ ...
    boundary;
    interiorPoints;
    axisPoint];

boundaryConstraints = [ ...
    (1:nBoundary).', ...
    [2:nBoundary,1].'];

meshP1 = delaunayTriangulation( ...
    P1points,boundaryConstraints);

triangleCenters = incenter(meshP1);

triangleInside = inpolygon( ...
    triangleCenters(:,1), ...
    triangleCenters(:,2), ...
    Rb,Zb);

if any(~triangleInside)
    error(['Some P1 triangles lie outside the analytic ', ...
           'Solov''ev boundary.']);
end

%% ============================================================
% 11. Upgrade P1 mesh to P3
% =============================================================

P3 = generate_P3_nodes(meshP1);

Ndof = size(P3.points,1);
Nt = size(P3.elements,1);

fprintf('\n');
fprintf('P1 vertices           = %d\n', ...
    size(P3.P1points,1));

fprintf('P1 triangles          = %d\n',Nt);

fprintf('P3 global DOFs        = %d\n',Ndof);

fprintf('P3 boundary DOFs      = %d\n', ...
    numel(P3.boundaryNodes));

assert(all(P3.P1points(:,1) > 0), ...
    'The verification mesh contains R <= 0.');

%% ============================================================
% 12. Precompute quadrature
% =============================================================

quadratureDegree = 10;

quad = precompute_quadrature( ...
    quadratureDegree);

fprintf('Quadrature degree     = %d\n', ...
    quadratureDegree);

fprintf('Quadrature points     = %d\n', ...
    numel(quad.w));

fprintf('Quadrature weight sum = %.16f\n', ...
    sum(quad.w));

%% ============================================================
% 13. Construct normalized profile arrays
% =============================================================

nw = 101;

profiles = solov.makeProfiles(nw);

dp_dpsiN = ...
    profiles.dp_dpsiN;

FdF_dpsiN = ...
    profiles.FdF_dpsiN;

%% ============================================================
% 14. Evaluate exact nodal solution
% =============================================================

Rnodes = P3.points(:,1);
Znodes = P3.points(:,2);

psiExact = ...
    solov.psi(Rnodes,Znodes);

psiNExact = ...
    solov.psiN(Rnodes,Znodes);

% This measures the mismatch between the curved analytic LCFS
% and the straight-edged affine mesh boundary.
geometryBoundaryMismatch = max(abs( ...
    psiExact(P3.boundaryNodes) ...
    - solov.psiBoundary));

fprintf('Polygon boundary mismatch = %.3e\n', ...
    geometryBoundaryMismatch);

fprintf('Exact nodal psiN range     = [%+.6e, %+.6e]\n', ...
    min(psiNExact),max(psiNExact));

%% ============================================================
% 15. Check the exact differential equation
% =============================================================

deltaStarExact = ...
    solov.deltaStarPsi(Rnodes,Znodes);

strongSourceExact = ...
    solov.strongSource(Rnodes,Znodes);

strongEquationError = max(abs( ...
    -deltaStarExact-strongSourceExact));

fprintf('Strong-equation error      = %.3e\n', ...
    strongEquationError);

assert(strongEquationError < 1.e-12);

%% ============================================================
% 16. Assemble stiffness matrix
% =============================================================

K = assemble_K_matrix(P3,quad,true);

globalSymmetryError = ...
    norm(K-K.','fro') ...
    /max(norm(K,'fro'),eps);

globalNullResidual = ...
    norm(K*ones(Ndof,1),inf) ...
    /max(norm(K,'fro'),eps);

fprintf('\n');
fprintf('K symmetry error      = %.3e\n', ...
    globalSymmetryError);

fprintf('K constant residual   = %.3e\n', ...
    globalNullResidual);

assert(globalSymmetryError < 1.e-12);
assert(globalNullResidual < 1.e-12);

%% ============================================================
% 17. Direct linear Solov'ev solve
%
% Since the Solov'ev profiles are constant, the source is
% linear once the exact dpsi is supplied.
% =============================================================

fExact = assemble_RHS( ...
    P3,psiExact,quad, ...
    dp_dpsiN,FdF_dpsiN, ...
    solov.dpsi);

[KII,rhsI,bc] = apply_BC( ...
    P3,K,fExact,solov.psiBoundary);

psiLinear = zeros(Ndof,1);

psiLinear(bc.boundaryNodes) = ...
    bc.boundaryValues;

psiLinear(bc.freeNodes) = ...
    KII\rhsI;

linearLhs = ...
    K(bc.freeNodes,:)*psiLinear;

linearResidual = ...
    linearLhs-fExact(bc.freeNodes);

linearResidualScale = max([ ...
    norm(linearLhs), ...
    norm(fExact(bc.freeNodes)), ...
    eps]);

linearRelativeResidual = ...
    norm(linearResidual)/linearResidualScale;

linearError = ...
    psiLinear-psiExact;

linearRelativeL2 = ...
    norm(linearError) ...
    /max(norm(psiExact),eps);

linearRelativeInfinity = ...
    norm(linearError,inf) ...
    /max(norm(psiExact,inf),eps);

linearBoundaryError = norm( ...
    psiLinear(bc.boundaryNodes) ...
    - bc.boundaryValues,inf);

fprintf('\n');
fprintf('========================================================\n');
fprintf(' Direct linear verification\n');
fprintf('========================================================\n');

fprintf('Relative residual      = %.3e\n', ...
    linearRelativeResidual);

fprintf('Nodal relative L2      = %.3e\n', ...
    linearRelativeL2);

fprintf('Relative infinity      = %.3e\n', ...
    linearRelativeInfinity);

fprintf('Boundary error         = %.3e\n', ...
    linearBoundaryError);

assert(linearRelativeResidual < 1.e-10);
assert(linearBoundaryError < 1.e-13);

%% ============================================================
% 18. Full Picard iteration
% =============================================================

initialAmplitudeFactor = 0.8;

psi0 = ...
    initialAmplitudeFactor*psiExact;

psi0(P3.boundaryNodes) = ...
    solov.psiBoundary;

options = struct();

options.omega = 0.5;
options.maxIterations = 100;

options.updateTolerance = 1.e-9;
options.residualTolerance = 1.e-9;

options.axisMode = solov.axisMode;

options.symmetryTolerance = 1.e-12;
options.fluxSpanTolerance = 1.e-12;

options.verbose = true;

[psiPicard,picardResult] = Picard_iteration( ...
    P3,quad, ...
    dp_dpsiN,FdF_dpsiN, ...
    psi0,options);

%% ============================================================
% 19. Nodal Picard errors
% =============================================================

picardError = ...
    psiPicard-psiExact;

picardRelativeNodalL2 = ...
    norm(picardError) ...
    /max(norm(psiExact),eps);

picardRelativeFreeL2 = ...
    norm(picardError(picardResult.freeNodes)) ...
    /max(norm(psiExact(picardResult.freeNodes)),eps);

picardRelativeInfinity = ...
    norm(picardError,inf) ...
    /max(norm(psiExact,inf),eps);

picardBoundaryError = norm( ...
    psiPicard(picardResult.boundaryNodes) ...
    - solov.psiBoundary,inf);

%% ============================================================
% 20. Magnetic-axis diagnostics
% =============================================================

freeNodes = ...
    picardResult.freeNodes;

switch solov.axisMode
    case 'max'
        [psiAxisNumerical,localAxisIndex] = ...
            max(psiPicard(freeNodes));

    case 'min'
        [psiAxisNumerical,localAxisIndex] = ...
            min(psiPicard(freeNodes));

    otherwise
        error('Unknown axis mode: %s',solov.axisMode);
end

axisNode = ...
    freeNodes(localAxisIndex);

axisPointNumerical = ...
    P3.points(axisNode,:);

relativeAxisFluxError = ...
    abs(psiAxisNumerical-solov.psiAxis) ...
    /max(abs(solov.psiAxis),eps);

axisPositionError = ...
    norm(axisPointNumerical-solov.axisPoint);

%% ============================================================
% 21. Integrated L2 and H1-seminorm errors
% =============================================================

L2ErrorSquared = 0.0;
L2ExactSquared = 0.0;

H1ErrorSquared = 0.0;
H1ExactSquared = 0.0;

L1 = 1-quad.xi-quad.eta;
L2 = quad.xi;
L3 = quad.eta;

for e = 1:Nt
    geomNodes = ...
        P3.P1elements(e,:);

    Ve = ...
        P3.P1points(geomNodes,:);

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

    aJ = J(1,1);
    bJ = J(1,2);
    cJ = J(2,1);
    dJ = J(2,2);

    invJT = (1/detJ)*[ ...
         dJ, -cJ;
        -bJ,  aJ];

    Rq = ...
        L1*R1+L2*R2+L3*R3;

    Zq = ...
        L1*Z1+L2*Z2+L3*Z3;

    ids = ...
        P3.elements(e,:);

    psiElement = ...
        psiPicard(ids);

    psiNumericalQ = ...
        quad.N*psiElement;

    dNdR = ...
        invJT(1,1)*quad.dNdxi ...
        + invJT(1,2)*quad.dNdeta;

    dNdZ = ...
        invJT(2,1)*quad.dNdxi ...
        + invJT(2,2)*quad.dNdeta;

    dpsiNumerical_dR = ...
        dNdR*psiElement;

    dpsiNumerical_dZ = ...
        dNdZ*psiElement;

    psiExactQ = ...
        solov.psi(Rq,Zq);

    dpsiExact_dR = ...
        solov.dpsi_dR(Rq,Zq);

    dpsiExact_dZ = ...
        solov.dpsi_dZ(Rq,Zq);

    measure = ...
        quad.w*detJ;

    fluxErrorQ = ...
        psiNumericalQ-psiExactQ;

    L2ErrorSquared = ...
        L2ErrorSquared ...
        + sum(measure.*abs(fluxErrorQ).^2);

    L2ExactSquared = ...
        L2ExactSquared ...
        + sum(measure.*abs(psiExactQ).^2);

    gradientErrorSquared = ...
          abs(dpsiNumerical_dR-dpsiExact_dR).^2 ...
        + abs(dpsiNumerical_dZ-dpsiExact_dZ).^2;

    exactGradientSquared = ...
          abs(dpsiExact_dR).^2 ...
        + abs(dpsiExact_dZ).^2;

    H1ErrorSquared = ...
        H1ErrorSquared ...
        + sum(measure.*gradientErrorSquared);

    H1ExactSquared = ...
        H1ExactSquared ...
        + sum(measure.*exactGradientSquared);
end

relativeL2Error = ...
    sqrt(L2ErrorSquared) ...
    /max(sqrt(L2ExactSquared),eps);

relativeH1Error = ...
    sqrt(H1ErrorSquared) ...
    /max(sqrt(H1ExactSquared),eps);

%% ============================================================
% 22. Final verification report
% =============================================================

fprintf('\n');
fprintf('========================================================\n');
fprintf(' General Solov''ev verification result\n');
fprintf('========================================================\n');

fprintf('Picard converged            = %d\n', ...
    picardResult.converged);

fprintf('Picard iterations           = %d\n', ...
    picardResult.iterations);

fprintf('Final nonlinear residual    = %.3e\n', ...
    picardResult.residualHistory(end));

fprintf('Final relative update       = %.3e\n', ...
    picardResult.updateHistory(end));

fprintf('Nodal relative L2 error     = %.3e\n', ...
    picardRelativeNodalL2);

fprintf('Free-node relative L2 error = %.3e\n', ...
    picardRelativeFreeL2);

fprintf('Relative infinity error     = %.3e\n', ...
    picardRelativeInfinity);

fprintf('Integrated relative L2      = %.3e\n', ...
    relativeL2Error);

fprintf('Integrated relative H1      = %.3e\n', ...
    relativeH1Error);

fprintf('Boundary error              = %.3e\n', ...
    picardBoundaryError);

fprintf('Exact psi_axis              = %+14.7e\n', ...
    solov.psiAxis);

fprintf('Numerical psi_axis          = %+14.7e\n', ...
    psiAxisNumerical);

fprintf('Relative axis-flux error    = %.3e\n', ...
    relativeAxisFluxError);

fprintf('Exact axis location         = (%.8f, %.8f)\n', ...
    solov.axisPoint(1),solov.axisPoint(2));

fprintf('Numerical axis node         = (%.8f, %.8f)\n', ...
    axisPointNumerical(1),axisPointNumerical(2));

fprintf('Axis-position error         = %.3e\n', ...
    axisPositionError);

fprintf('Polygon boundary mismatch   = %.3e\n', ...
    geometryBoundaryMismatch);

fprintf('========================================================\n');

assert(picardResult.converged, ...
    'Picard iteration did not converge.');

assert(picardBoundaryError < 1.e-13, ...
    'Picard solution violates the fixed boundary condition.');

assert(picardResult.residualHistory(end) ...
    < options.residualTolerance, ...
    'Final Picard residual is too large.');

%% ============================================================
% 23. Construct P3 visualization connectivity
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

%% ============================================================
% 24. Plot LCFS and mesh
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
    solov.axisPoint(1), ...
    solov.axisPoint(2), ...
    'ro', ...
    'MarkerFaceColor','r');

plot( ...
    [Rinner,Router,Rtop], ...
    [0,0,Ztop], ...
    'bs', ...
    'MarkerFaceColor','b');

xlabel('R');
ylabel('Z');

legend( ...
    'P1 mesh', ...
    'analytic LCFS', ...
    'magnetic axis', ...
    'shape control points', ...
    'Location','best');

title('General Solov''ev LCFS and mesh');

%% ============================================================
% 25. Plot exact, numerical, and error fields
% =============================================================

figure('Color','w');

subplot(1,3,1);

trisurf( ...
    Tplot, ...
    Rnodes,Znodes,psiExact, ...
    'EdgeColor','none');

view(2);
colorbar;

xlabel('R');
ylabel('Z');
title('Exact Solov''ev \psi');

subplot(1,3,2);

trisurf( ...
    Tplot, ...
    Rnodes,Znodes,psiPicard, ...
    'EdgeColor','none');

view(2);
colorbar;

xlabel('R');
ylabel('Z');
title('P3 Picard solution');

subplot(1,3,3);

trisurf( ...
    Tplot, ...
    Rnodes,Znodes,picardError, ...
    'EdgeColor','none');

view(2);
colorbar;

xlabel('R');
ylabel('Z');
title('\psi_h-\psi_{exact}');

sgtitle('General P3 Solov''ev verification');

%% ============================================================
% 26. Plot Picard convergence
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

title('General Solov''ev Picard convergence');