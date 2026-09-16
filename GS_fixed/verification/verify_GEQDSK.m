%% verify_GEQDSK.m
% Cross-validate solve_GS against one GEQDSK equilibrium.
clearvars; close all; clc;

thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir,'..','..');
addpath(thisDir,fullfile(rootDir,'equilibrium'));

%% Reference equilibrium and FEM controls
gfile = fullfile(rootDir,'gfiles','geqdsk_PT0.6');
NR = 128; NZ = 128;
quadratureDegree = 10;
nCompareR = 181;

reference = read_geqdsk(gfile);
eqfunc = build_interpolants(reference);
boundary = clean_boundary(reference.rbbbs,reference.zbbbs);

psiSpan = reference.sibry-reference.simag;
if ~isfinite(psiSpan) || psiSpan==0
    error('VERIFY_GEQDSK:FluxSpan','GEQDSK flux span must be nonzero.');
end

% solve_GS uses psi=0 on the boundary and a positive axis value.  The
% transformed GEQDSK reference has the same psiN and this native convention.
fluxSign = sign(psiSpan);
pprimeN = reference.pprime(:)*psiSpan;
FFprimeN = reference.ffprim(:)*psiSpan;

%% Solve with the GEQDSK boundary and profiles
input = struct();
input.grid = struct('nBoundary',size(boundary,1),'NR',NR,'NZ',NZ, ...
    'quadratureDegree',quadratureDegree,'nProfile',reference.nw);
input.dim = struct('R0',reference.rcentr, ...
    'B0',reference.fpol(end)/reference.rcentr, ...
    'aMinor',(max(boundary(:,1))-min(boundary(:,1)))/2);
input.profile = struct('mode',1,'pprime',pprimeN, ...
    'FFprime',FFprimeN,'pedge',reference.pres(end));
input.boundary = struct('R',boundary(:,1),'Z',boundary(:,2), ...
    'psiBoundary',0);
input.Picard = struct('omega',0.5,'maxIterations',150, ...
    'updateTolerance',1e-9,'residualTolerance',1e-9, ...
    'axisMode','max','symmetryTolerance',1e-12, ...
    'fluxSpanTolerance',1e-12);
input.output = struct('verbose',true,'checkTime',true,'showPlot',false);

equilibrium = solve_GS(input);

%% Compare psi on a common uniform grid inside the FEM domain
Rvec = linspace(min(boundary(:,1)),max(boundary(:,1)),nCompareR);
aspect = range(boundary(:,2))/range(boundary(:,1));
Zvec = linspace(min(boundary(:,2)),max(boundary(:,2)), ...
    max(121,round(nCompareR*aspect)));
[R,Z] = meshgrid(Rvec,Zvec);

P3 = equilibrium.mesh;
P1 = triangulation(P3.P1elements,P3.P1points);
psiFEM = reshape(eval_P3_sol(P1,P3,equilibrium.psi,[R(:),Z(:)]),size(R));
psiGEQDSK = reshape(eqfunc.psi(R(:),Z(:)),size(R));
psiReference = fluxSign*(reference.sibry-psiGEQDSK);

valid = isfinite(psiFEM) & isfinite(psiReference);
errorPsi = psiFEM(valid)-psiReference(valid);
referencePsi = psiReference(valid);
psiNFEM = (psiFEM-equilibrium.solverResult.psiAxis) ...
    /equilibrium.solverResult.dpsi;
psiNGEQDSK = (psiGEQDSK-reference.simag)/psiSpan;
errorPsiN = psiNFEM(valid)-psiNGEQDSK(valid);

relativeL2 = norm(errorPsi)/max(norm(referencePsi),eps);
relativeInfinity = norm(errorPsi,inf)/max(norm(referencePsi,inf),eps);
normalizedRMS = sqrt(mean((errorPsi/abs(psiSpan)).^2));
psiNRMS = sqrt(mean(errorPsiN.^2));
psiNInfinity = norm(errorPsiN,inf);
axisFluxError = abs(equilibrium.solverResult.psiAxis-abs(psiSpan))/abs(psiSpan);
axisPositionError = norm(equilibrium.axisPoint- ...
    [reference.rmaxis,reference.zmaxis]);

referenceBoundaryPsi = eqfunc.psi(boundary(:,1),boundary(:,2));
geometryBoundaryMismatch = norm( ...
    referenceBoundaryPsi-reference.sibry,inf)/abs(psiSpan);
boundaryConditionError = norm(equilibrium.psi(P3.boundaryNodes),inf);

verification = struct('gfile',gfile,'relativeL2',relativeL2, ...
    'relativeInfinity',relativeInfinity,'normalizedRMS',normalizedRMS, ...
    'psiNRMS',psiNRMS,'psiNInfinity',psiNInfinity, ...
    'axisFluxRelativeError',axisFluxError, ...
    'axisPositionError',axisPositionError, ...
    'geometryBoundaryMismatch',geometryBoundaryMismatch, ...
    'boundaryConditionError',boundaryConditionError, ...
    'nComparisonPoints',nnz(valid));

fprintf('\nGEQDSK cross-validation: %s\n',gfile);
fprintf('  comparison points             = %d\n',nnz(valid));
fprintf('  relative L2 psi error         = %.3e\n',relativeL2);
fprintf('  relative infinity psi error   = %.3e\n',relativeInfinity);
fprintf('  RMS error / |Delta psi|       = %.3e\n',normalizedRMS);
fprintf('  psiN RMS / infinity error     = %.3e / %.3e\n', ...
    psiNRMS,psiNInfinity);
fprintf('  relative axis-flux error      = %.3e\n',axisFluxError);
fprintf('  axis-position error [m]       = %.3e\n',axisPositionError);
fprintf('  GEQDSK polygon mismatch       = %.3e\n',geometryBoundaryMismatch);
fprintf('  FEM boundary-condition error  = %.3e\n',boundaryConditionError);

assert(equilibrium.solverResult.converged, ...
    'VERIFY_GEQDSK:Convergence','Picard iteration did not converge.');
assert(all(isfinite(errorPsi)) && nnz(valid)>100, ...
    'VERIFY_GEQDSK:Comparison','Insufficient finite comparison points.');
assert(boundaryConditionError<1e-13, ...
    'VERIFY_GEQDSK:BoundaryCondition','FEM boundary condition is violated.');
assert(relativeL2<2e-3 && relativeInfinity<5e-3, ...
    'VERIFY_GEQDSK:FluxError','FEM and GEQDSK psi fields do not agree.');
assert(axisFluxError<5e-3 && axisPositionError<2e-2, ...
    'VERIFY_GEQDSK:AxisError','FEM and GEQDSK magnetic axes do not agree.');
assert(geometryBoundaryMismatch<1e-3, ...
    'VERIFY_GEQDSK:BoundaryGeometry', ...
    'The polygonal boundary is not a GEQDSK flux surface.');
assert(isequal(equilibrium.profile.dp_dpsiN,pprimeN) && ...
    isequal(equilibrium.profile.FdF_dpsiN,FFprimeN), ...
    'VERIFY_GEQDSK:Profiles','Normalized GEQDSK profiles were not preserved.');

plot_comparison(R,Z,psiFEM,psiGEQDSK,errorPsi,valid, ...
    equilibrium,reference,boundary,psiSpan);

%% Local functions
function boundary = clean_boundary(R,Z)
boundary = [R(:),Z(:)];
if size(boundary,1)<3 || any(~isfinite(boundary),'all')
    error('VERIFY_GEQDSK:Boundary','GEQDSK boundary is invalid.');
end
scale = max([range(boundary(:,1)),range(boundary(:,2)),1]);
tol = 1e-12*scale;
keep = [true; hypot(diff(boundary(:,1)),diff(boundary(:,2)))>tol];
boundary = boundary(keep,:);
if norm(boundary(end,:)-boundary(1,:))<=tol, boundary(end,:) = []; end
if size(boundary,1)<3
    error('VERIFY_GEQDSK:Boundary','GEQDSK boundary has too few unique points.');
end
end

function plot_comparison(R,Z,psiFEM,psiGEQDSK,errorPsi,valid, ...
    equilibrium,reference,boundary,psiSpan)
psiNFEM = (psiFEM-equilibrium.solverResult.psiAxis) ...
    /equilibrium.solverResult.dpsi;
psiNGEQDSK = (psiGEQDSK-reference.simag)/psiSpan;
psiNGEQDSK(~valid) = NaN;
errorField = nan(size(R));
errorField(valid) = errorPsi/abs(psiSpan);

figure('Color','w','Name','P3 GS and GEQDSK comparison');
layout = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax(1) = nexttile(layout);
contourf(ax(1),R,Z,psiNFEM,0:.05:1,'LineStyle','none');
title(ax(1),'FEM \psi_N'); clim(ax(1),[0,1]); colorbar(ax(1));
ax(2) = nexttile(layout);
contourf(ax(2),R,Z,psiNGEQDSK,0:.05:1,'LineStyle','none');
title(ax(2),'GEQDSK \psi_N'); clim(ax(2),[0,1]); colorbar(ax(2));
ax(3) = nexttile(layout);
contourf(ax(3),R,Z,errorField,31,'LineStyle','none');
limit = max(abs(errorField),[],'all','omitnan');
if limit>0, clim(ax(3),[-limit,limit]); end
title(ax(3),'(\psi_{FEM}-\psi_{GEQDSK})/|\Delta\psi|'); colorbar(ax(3));
for a = ax
    hold(a,'on'); plot(a,boundary(:,1),boundary(:,2),'k-','LineWidth',1.2);
    axis(a,'equal'); grid(a,'on'); box(a,'on'); xlabel(a,'R [m]'); ylabel(a,'Z [m]');
end

result = equilibrium.solverResult;
figure('Color','w','Name','GEQDSK cross-validation convergence');
iteration = 1:result.iterations;
semilogy(iteration,result.updateHistory,'o-',iteration, ...
    result.residualHistory,'s-','LineWidth',1.2);
grid on; box on; xlabel('Picard iteration'); ylabel('Relative error');
legend('solution update','nonlinear residual','Location','best');
end
