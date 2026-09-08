%% Fixed-boundary GS reconstruction of geqdsk_PT0.6
% The LCFS, p', FF', and F profiles are imported from the reference
% GEQDSK. Only files inside GS_fixed are modified; the existing sibling
% GEQDSK reader is called as a dependency.

clear
clc
close all

projectRoot = fileparts(mfilename('fullpath'));
parentRoot = fileparts(projectRoot);
addpath(fullfile(projectRoot,'geometry'), ...
        fullfile(projectRoot,'fem'), ...
        fullfile(projectRoot,'physics'), ...
        fullfile(projectRoot,'solver'), ...
        fullfile(projectRoot,'post'), ...
        fullfile(projectRoot,'ballooning'), ...
        fullfile(parentRoot,'equilibrium'));

%% Reference GEQDSK and user configuration
geqdskFile = fullfile(parentRoot,'Miller','geqdsk_PT0.6');
if ~isfile(geqdskFile)
    error('GS:mainPT06:MissingGEQDSK', ...
        'Reference file was not found: %s',geqdskFile);
end
reference = read_geqdsk(geqdskFile);
psiNProfile = linspace(0,1,reference.nw).';

boundaryR = reference.rbbbs(:);
boundaryZ = reference.zbbbs(:);
Rout = max(boundaryR);
Rin = min(boundaryR);
R0 = 0.5*(Rout+Rin);
minorRadius = 0.5*(Rout-Rin);
verticalHalfHeight = 0.5*(max(boundaryZ)-min(boundaryZ));

% Map each ordered GEQDSK point to a corresponding point on a circle. This
% supplies a safe fixed-boundary homotopy without replacing the final LCFS
% by a Miller approximation.
geometricAngle = unwrap(atan2( ...
    boundaryZ/verticalHalfHeight,(boundaryR-R0)/minorRadius));
startR = R0+minorRadius*cos(geometricAngle);
startZ = minorRadius*sin(geometricAngle);

cfg = struct();
cfg.referenceFile = geqdskFile;
boundaryCount = numel(boundaryR);
if hypot(boundaryR(end)-boundaryR(1),boundaryZ(end)-boundaryZ(1))<1e-12
    boundaryCount = boundaryCount-1;
end
cfg.startBoundary = struct('type','tabulated', ...
    'R',startR,'Z',startZ,'nBoundary',boundaryCount);
cfg.targetBoundary = struct('type','tabulated', ...
    'R',boundaryR,'Z',boundaryZ,'nBoundary',boundaryCount);
cfg.continuationLambdas = (0:0.2:1).';
cfg.mesh = struct('targetH',0.02,'boundaryClearance',0.12);

pSpec = struct('type','tabulated','psiN',psiNProfile, ...
    'values',reference.pprime(:),'method','pchip');
FFspec = struct('type','tabulated','psiN',psiNProfile, ...
    'values',reference.ffprim(:),'method','pchip', ...
    'FValues',reference.fpol(:));
cfg.profiles = make_profiles(struct( ...
    'flux',struct('psiAxis',reference.simag, ...
                  'psiBoundary',reference.sibry, ...
                  'rangePolicy','error','tolerance',1e-8), ...
    'pprime',pSpec,'FFprime',FFspec));

cfg.initialGuess = struct('boundaryValue',reference.sibry, ...
    'axisValue',reference.simag);
cfg.picard = struct( ...
    'tolerance',1e-8, ...
    'residualTolerance',2e-8, ...
    'maxIterations',2000, ...
    'omega',0.15, ...
    'boundaryValue',reference.sibry, ...
    'axisMode','min', ...
    'updateFluxAxis',true, ...
    'failOnNonconvergence',true, ...
    'verbose',false);
cfg.topologyLevels = (0.1:0.1:0.9).';
cfg.FBoundary = reference.fpol(end);
cfg.makePlots = true;

%% Nonlinear GS solve
continuationOptions = struct('mesh',cfg.mesh, ...
    'initialGuess',cfg.initialGuess,'picard',cfg.picard,'verbose',true);
continuation = solve_continuation( ...
    cfg.startBoundary,cfg.targetBoundary,cfg.continuationLambdas, ...
    cfg.profiles,continuationOptions);

geom = continuation.finalGeometry;
mesh = continuation.finalMesh;
psi = continuation.finalPsi;
solverResult = continuation.finalSolver;
axisData = find_axis(mesh,psi,struct('mode','min', ...
    'boundaryValue',cfg.picard.boundaryValue));
topology = check_topology(mesh,psi,axisData,struct( ...
    'normalizedLevels',cfg.topologyLevels, ...
    'boundaryValue',cfg.picard.boundaryValue, ...
    'failOnFailure',false,'requireSingleO',true,'requireNoX',true));

finalFluxSpec = cfg.profiles.flux;
finalFluxSpec.psiAxis = axisData.psi;
finalFluxSpec.psiBoundary = cfg.picard.boundaryValue;
FofPsi = @(psiValue) reconstruct_F( ...
    psiValue,cfg.FBoundary,cfg.profiles.FFprime,finalFluxSpec);
field = compute_B(mesh,psi,FofPsi);

equilibrium = struct();
equilibrium.configuration = cfg;
equilibrium.geometry = geom;
equilibrium.mesh = mesh;
equilibrium.psi = psi;
equilibrium.profiles = cfg.profiles;
equilibrium.solver = solverResult;
equilibrium.continuation = continuation;
equilibrium.axis = axisData;
equilibrium.topology = topology;
equilibrium.field = field;
equilibrium.referenceGEQDSK = reference;

%% Interior q-profile and integral-current validation
validationPsiN = (0.1:0.1:0.9).';
bundle = extract_local_surface_bundle(mesh,psi,axisData,validationPsiN, ...
    struct('nPoloidal',512,'boundaryValue',cfg.picard.boundaryValue));
localEq = build_local_pest_equilibrium( ...
    equilibrium,bundle,struct('nTheta',512));
equilibrium.localEquilibrium = localEq;

qReference = interp1(psiNProfile,reference.qpsi(:), ...
    validationPsiN,'pchip');
qFEM = localEq.profiles.q;
qRelativeRms = norm(qFEM-qReference)/norm(qReference);
qComparison = table(validationPsiN,qReference,qFEM, ...
    (qFEM-qReference)./qReference, ...
    'VariableNames',{'psiN','qGEQDSK','qFEM','relativeError'});

plasmaCurrent = integrate_plasma_current( ...
    mesh,psi,cfg.profiles,finalFluxSpec);
currentMagnitudeError = (abs(plasmaCurrent)-abs(reference.current)) ...
    /abs(reference.current);

shape = measure_boundary_shape(geom.R,geom.Z);
referenceShape = measure_boundary_shape(boundaryR,boundaryZ);
axisPositionError = hypot(axisData.R-reference.rmaxis, ...
    axisData.Z-reference.zmaxis);
axisFluxRelativeError = (axisData.psi-reference.simag) ...
    /abs(reference.simag);

validation = struct('qComparison',qComparison, ...
    'qRelativeRms',qRelativeRms,'plasmaCurrent',plasmaCurrent, ...
    'currentMagnitudeRelativeError',currentMagnitudeError, ...
    'shape',shape,'referenceShape',referenceShape, ...
    'axisPositionError',axisPositionError, ...
    'axisFluxRelativeError',axisFluxRelativeError);
equilibrium.validation = validation;

%% Text summary
fprintf('\ngeqdsk_PT0.6 fixed-boundary reconstruction complete.\n')
fprintf('Reference file       : %s\n',geqdskFile)
fprintf('Final mesh           : %d nodes, %d P1 elements\n', ...
    mesh.nNodes,mesh.nElements)
fprintf('Picard solve         : %d iterations, residual = %.3e\n', ...
    solverResult.nIterations,solverResult.finalResidual)
fprintf('Topology             : passed=%d, O=%d, X=%d\n', ...
    topology.passed,topology.nOPoints,topology.nXPoints)
fprintf(['Magnetic axis        : FEM=(%.8f, %.8f), ' ...
         'GEQDSK=(%.8f, %.8f), |Delta|=%.3e m\n'], ...
    axisData.R,axisData.Z,reference.rmaxis,reference.zmaxis, ...
    axisPositionError)
fprintf(['Axis flux            : FEM=%.9g, GEQDSK=%.9g, ' ...
         'relative error=%+.3e\n'], ...
    axisData.psi,reference.simag,axisFluxRelativeError)
fprintf(['Plasma current       : |FEM|=%.6f MA, GEQDSK=%.6f MA, ' ...
         'relative error=%+.3e\n'], ...
    abs(plasmaCurrent)/1e6,abs(reference.current)/1e6,currentMagnitudeError)
fprintf('q-profile relative RMS error = %.3e\n',qRelativeRms)
fprintf(['LCFS shape           : R0=%.6f m, a=%.6f m, ' ...
         'kappa=%.6f, delta=%.6f\n\n'], ...
    shape.R0,shape.a,shape.kappa,shape.delta)
disp(qComparison)

%% Plots
if cfg.makePlots
    plot_flux(mesh,psi,struct('axisData',axisData, ...
        'normalizedLevels',cfg.topologyLevels, ...
        'boundaryValue',cfg.picard.boundaryValue, ...
        'title','GEQDSK PT0.6 fixed-boundary reconstruction'));

    figure('Name','PT0.6 boundary comparison','Color','w')
    plot(reference.rbbbs,reference.zbbbs,'k-','LineWidth',2, ...
        'DisplayName','GEQDSK LCFS')
    hold on
    plot(geom.closedPoints(:,1),geom.closedPoints(:,2),'r--', ...
        'LineWidth',1.4,'DisplayName','GS fixed LCFS')
    hold off
    axis equal tight
    grid on
    xlabel('R [m]')
    ylabel('Z [m]')
    legend('Location','best')

    figure('Name','PT0.6 q comparison','Color','w')
    plot(psiNProfile,reference.qpsi,'k-','LineWidth',2, ...
        'DisplayName','GEQDSK')
    hold on
    plot(validationPsiN,qFEM,'ro-','LineWidth',1.4, ...
        'DisplayName','GS fixed FEM')
    hold off
    grid on
    xlabel('\psi_N')
    ylabel('q')
    legend('Location','best')

    figure('Name','PT0.6 source profiles','Color','w')
    tiledlayout(1,2)
    nexttile
    plot(psiNProfile,reference.pprime,'LineWidth',1.8)
    grid on
    xlabel('\psi_N')
    ylabel('p'' [Pa/(Wb/rad)]')
    nexttile
    plot(psiNProfile,reference.ffprim,'LineWidth',1.8)
    grid on
    xlabel('\psi_N')
    ylabel('FF'' [T^2 m^2/(Wb/rad)]')
end


function current = integrate_plasma_current(mesh,psi,profiles,fluxSpec)
elements = mesh.elements;
nodes = mesh.nodes;
p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
area = 0.5*abs((p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) ...
    -(p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1)));
centroidR = (p1(:,1)+p2(:,1)+p3(:,1))/3;
psiElement = mean(reshape(psi(elements(:)),size(elements)),2);
pressureDerivative = pprime(psiElement,profiles.pprime,fluxSpec);
FFDerivative = FFprime(psiElement,profiles.FFprime,fluxSpec);
jphi = centroidR.*pressureDerivative ...
    +FFDerivative./(profiles.mu0*centroidR);
current = sum(jphi.*area);
end


function shape = measure_boundary_shape(R,Z)
R = R(:);
Z = Z(:);
[Rout,~] = max(R);
[Rin,~] = min(R);
R0 = 0.5*(Rout+Rin);
a = 0.5*(Rout-Rin);
[Ztop,topIndex] = max(Z);
[Zbottom,bottomIndex] = min(Z);
RtopBottom = 0.5*(R(topIndex)+R(bottomIndex));
shape = struct('R0',R0,'a',a,'aspectRatio',R0/a, ...
    'kappa',(Ztop-Zbottom)/(2*a), ...
    'delta',(R0-RtopBottom)/a, ...
    'Rout',Rout,'Rin',Rin,'Rtop',R(topIndex),'Ztop',Ztop);
end
