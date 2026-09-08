%% Fixed-boundary Grad-Shafranov solver for a concave C-shaped LCFS
% This script connects geometry, P1 FEM, nonlinear profiles, Picard
% continuation, and topology/magnetic-field post-processing.
%
% The default profile coefficients form a robust nonlinear reference case;
% they are not a fitted equilibrium for a particular experiment. Replace
% cfg.profiles with experimentally constrained p'(psi) and FF'(psi) before
% interpreting pressure, current, q, or absolute field values physically.

clear
clc
close all

projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot,'geometry'), ...
        fullfile(projectRoot,'fem'), ...
        fullfile(projectRoot,'physics'), ...
        fullfile(projectRoot,'solver'), ...
        fullfile(projectRoot,'post'));

%% User configuration
cfg = struct();

% Pointed boomerang indentation. Two power-law arcs meet at
% (Rend,+/-Zend); at Z=0 their crossings are Rback and Rtip. The target
% has Rtip<A while both pointed ends sweep outward to Rend.
cfg.startBoundary = struct( ...
    'type','boomerang-c','A',3.0, ...
    'Rback',2.0,'Rtip',3.4,'Rend',3.8,'Zend',1.3, ...
    'pInner',1.6,'pOuter',1.4,'nBoundary',200);
cfg.targetBoundary = struct( ...
    'type','boomerang-c','A',3.0, ...
    'Rback',2.0,'Rtip',2.65,'Rend',3.8,'Zend',2.7, ...
    'pInner',1.6,'pOuter',1.4,'nBoundary',400);
cfg.continuationLambdas = (0:0.05:1).';

% Coarse grid
%cfg.mesh = struct('targetH',0.11,'boundaryClearance',0.12);
% Intermediate grid
%cfg.mesh = struct('targetH',0.08,'boundaryClearance',0.12);
% Dense grid
cfg.mesh = struct('targetH',0.05,'boundaryClearance',0.12);

% Negative source with psi_LCFS=0 gives the convention psi_axis<0.
cfg.profiles = make_profiles(struct( ...
    'pprime',struct('type','normalized-power', ...
                    'axisValue',-8.0e5,'exponent',2.0), ...
    'FFprime',struct('type','normalized-power', ...
                     'axisValue',-0.8,'exponent',1.0)));

cfg.initialGuess = struct('boundaryValue',0.0,'axisValue',-1.0);
cfg.picard = struct( ...
    'tolerance',1e-8, ...
    'residualTolerance',2e-8, ...
    'maxIterations',600, ...
    'omega',0.3, ...
    'boundaryValue',0.0, ...
    'axisMode','min', ...
    'updateFluxAxis',true, ...
    'failOnNonconvergence',true, ...
    'verbose',false);

cfg.topologyLevels = (0.1:0.1:0.9).';
cfg.FBoundary = 3.0;  % F=R*B_phi at the LCFS for field reconstruction
cfg.makePlots = true;

%% Follow the equilibrium branch into the target C-shape
continuationOptions = struct( ...
    'mesh',cfg.mesh, ...
    'initialGuess',cfg.initialGuess, ...
    'picard',cfg.picard, ...
    'verbose',true);

continuation = solve_continuation( ...
    cfg.startBoundary,cfg.targetBoundary, ...
    cfg.continuationLambdas,cfg.profiles,continuationOptions);

geom = continuation.finalGeometry;
mesh = continuation.finalMesh;
psi = continuation.finalPsi;
solverResult = continuation.finalSolver;

%% Axis and topology diagnostics
axisData = find_axis(mesh,psi,struct( ...
    'mode',cfg.picard.axisMode, ...
    'boundaryValue',cfg.picard.boundaryValue));

topology = check_topology(mesh,psi,axisData,struct( ...
    'normalizedLevels',cfg.topologyLevels, ...
    'boundaryValue',cfg.picard.boundaryValue, ...
    'failOnFailure',false, ...
    'requireSingleO',true, ...
    'requireNoX',true));

if ~topology.passed
    warning('GS:main:TopologyCheckFailed', ...
        ['The FEM solve converged, but the requested single-axis topology ' ...
         'check failed. Inspect equilibrium.topology before using the result.']);
end

%% Magnetic field
% The normalized-power FF' profile can be integrated analytically from the
% LCFS value of F. This keeps the reconstructed B_phi consistent with the
% FF' used in the Grad-Shafranov source.
FFspec = cfg.profiles.FFprime;
FofPsi = @(psiValue) reconstruct_F(psiValue,axisData.psi, ...
    cfg.picard.boundaryValue,cfg.FBoundary, ...
    FFspec.axisValue,FFspec.exponent);
field = compute_B(mesh,psi,FofPsi);
indentation = struct();
indentation.tipR = geom.outboardMidplaneR;
indentation.nominalMajorRadius = geom.nominalMajorRadius;
indentation.magneticAxisR = axisData.R;
indentation.tipInsideNominalMajorRadius = ...
    geom.outboardMidplaneR < geom.nominalMajorRadius;
indentation.tipInsideMagneticAxis = geom.outboardMidplaneR < axisData.R;

%% Collect the complete result in one workspace variable
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
equilibrium.indentation = indentation;

%% Text summary
nSteps = numel(continuation.steps);
lambda = zeros(nSteps,1);
nNodes = zeros(nSteps,1);
nElements = zeros(nSteps,1);
nPicard = zeros(nSteps,1);
finalResidual = zeros(nSteps,1);
psiAxis = zeros(nSteps,1);
for k = 1:nSteps
    lambda(k) = continuation.steps(k).lambda;
    nNodes(k) = continuation.steps(k).mesh.nNodes;
    nElements(k) = continuation.steps(k).mesh.nElements;
    nPicard(k) = continuation.steps(k).solver.nIterations;
    finalResidual(k) = continuation.steps(k).solver.finalResidual;
    psiAxis(k) = continuation.steps(k).solver.axisValue;
end
continuationTable = table(lambda,nNodes,nElements,nPicard, ...
    finalResidual,psiAxis);

fprintf('\nFixed-boundary Grad-Shafranov solve complete.\n')
disp(continuationTable)
fprintf('Final magnetic axis : R = %.8f, Z = %.8f, psi = %.8g\n', ...
    axisData.R,axisData.Z,axisData.psi)
fprintf('Final Picard solve  : %d iterations, residual = %.3e\n', ...
    solverResult.nIterations,solverResult.finalResidual)
fprintf('Flux topology       : passed = %d, O-points = %d, X-points = %d\n', ...
    topology.passed,topology.nOPoints,topology.nXPoints)
fprintf('Final mesh          : %d nodes, %d P1 elements\n\n', ...
    mesh.nNodes,mesh.nElements)
fprintf(['Boomerang geometry  : R_tip = %.8f, R_end = %.8f, ' ...
         'nominal R0 = %.8f, elongation = %.6f\n\n'], ...
    geom.outboardMidplaneR,geom.endPointR,geom.nominalMajorRadius, ...
    geom.boundingBoxElongation)
fprintf(['Indentation checks  : tip < nominal R0 = %d, ' ...
         'tip < magnetic-axis R = %d\n\n'], ...
    indentation.tipInsideNominalMajorRadius, ...
    indentation.tipInsideMagneticAxis)

%% Plots
if cfg.makePlots
    plot_flux(mesh,psi,struct( ...
        'axisData',axisData, ...
        'normalizedLevels',cfg.topologyLevels, ...
        'boundaryValue',cfg.picard.boundaryValue, ...
        'title','Final pointed boomerang C-shape equilibrium'));

    figure('Name','Picard convergence','Color','w')
    semilogy(1:solverResult.nIterations,solverResult.updateHistory, ...
        '-','LineWidth',1.6)
    hold on
    semilogy(1:solverResult.nIterations,solverResult.residualHistory, ...
        '--','LineWidth',1.6)
    hold off
    grid on
    xlabel('Picard iteration')
    ylabel('Relative measure')
    legend('iterate update','equation residual','Location','best')
    title('Final continuation-step convergence')

    figure('Name','C-shape continuation','Color','w')
    colorMap = parula(nSteps);
    hold on
    for k = 1:nSteps
        boundary = continuation.steps(k).geom.closedPoints;
        plot(boundary(:,1),boundary(:,2),'-', ...
            'Color',colorMap(k,:),'LineWidth',1.4, ...
            'DisplayName',sprintf('lambda = %.1f',lambda(k)))
    end
    hold off
    axis equal tight
    grid on
    xlabel('R')
    ylabel('Z')
    title('LCFS continuation into pointed boomerang shaping')
    legend('Location','bestoutside')

    figure('Name','Poloidal magnetic field','Color','w')
    trisurf(mesh.elements,mesh.nodes(:,1),mesh.nodes(:,2), ...
        field.node.Bp,'EdgeColor','none','FaceColor','interp')
    view(2)
    axis equal tight
    grid on
    colorbar
    xlabel('R')
    ylabel('Z')
    title('|B_p| on the final P1 mesh')
end


function F = reconstruct_F(psi,psiAxis,psiBoundary,FBoundary, ...
        FFaxisValue,exponent)
% Integrate FF'=F*dF/dpsi for FF'=C*(1-psi_N)^exponent.
span = psiBoundary-psiAxis;
psiNormalized = (psi-psiAxis)/span;
psiNormalized = min(max(psiNormalized,0),1);
F2 = FBoundary^2 ...
   -2*FFaxisValue*span/(exponent+1) ...
    .*(1-psiNormalized).^(exponent+1);
if any(F2(:)<=0)
    error('GS:main:NonPositiveFSquared', ...
        ['The selected FBoundary and FFprime profile produce F^2<=0. ' ...
         'Increase FBoundary or reduce the FFprime amplitude.']);
end
F = sqrt(F2);
end
