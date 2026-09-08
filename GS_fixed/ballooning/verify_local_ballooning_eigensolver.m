function report = verify_local_ballooning_eigensolver(input,doPlot)
%VERIFY_LOCAL_BALLOONING_EIGENSOLVER Verify sparse local mode convergence.
%
%   REPORT = VERIFY_LOCAL_BALLOONING_EIGENSOLVER(EQUILIBRIUM,TRUE) checks
%   sparse-vs-dense agreement, angular resolution, field-line-domain
%   convergence, residuals, localization, and Dirichlet boundaries.

narginchk(1,2);
if nargin<2
    doPlot = true;
end
if isstruct(input) && all(isfield(input,{'map','eqfunc','bundle'}))
    localEq = input;
else
    equilibrium = input;
    bundle = extract_local_surface_bundle( ...
        equilibrium.mesh,equilibrium.psi,equilibrium.axis, ...
        (0.55:0.05:0.95).', ...
        struct('nPoloidal',512,'boundaryValue',0));
    localEq = build_local_pest_equilibrium( ...
        equilibrium,bundle,struct('nTheta',512));
end

targetPsiN = 0.75;
base = struct('theta0',0,'numEigenvalues',4, ...
    'poloidalModeCutoff',24,'eigenTolerance',1e-10, ...
    'maxIterations',3000);

% Cross-check the sparse generalized solve against the existing dense
% solver on the same deliberately small assembled system.
smallOptions = base;
smallOptions.nPeriods = 1;
smallOptions.elementsPerPeriod = 64;
small = solve_local_ballooning_mode(localEq,targetPsiN,smallOptions);
add_sturm_liouville_functions();
[lambdaDense,~] = calculate_eigenvalue( ...
    small.thetaDof,small.matrices.G,small.matrices.C, ...
    small.matrices.F,'dirichlet');
denseSparseError = abs(lambdaDense(1)-small.lambda(1)) ...
    /max(1,abs(lambdaDense(1)));

domainPeriods = [4;6;8;10];
domainLambda = zeros(size(domainPeriods));
domainResidual = zeros(size(domainPeriods));
domainEdgeAmplitude = zeros(size(domainPeriods));
for k = 1:numel(domainPeriods)
    runOptions = base;
    runOptions.nPeriods = domainPeriods(k);
    runOptions.elementsPerPeriod = 128;
    run = solve_local_ballooning_mode( ...
        localEq,targetPsiN,runOptions);
    domainLambda(k) = run.lambda(1);
    domainResidual(k) = run.diagnostics.leadingRelativeResidual;
    domainEdgeAmplitude(k) = ...
        run.diagnostics.leadingEdgeAmplitudeRatio;
end
domainTable = table(domainPeriods,domainLambda,domainResidual, ...
    domainEdgeAmplitude,'VariableNames', ...
    {'nPeriods','lambda1','relativeResidual','edgeAmplitudeRatio'});

elementsPerPeriod = [64;128;256];
resolutionLambda = zeros(size(elementsPerPeriod));
resolutionResidual = zeros(size(elementsPerPeriod));
reference = [];
for k = 1:numel(elementsPerPeriod)
    runOptions = base;
    runOptions.nPeriods = 10;
    runOptions.elementsPerPeriod = elementsPerPeriod(k);
    run = solve_local_ballooning_mode( ...
        localEq,targetPsiN,runOptions);
    resolutionLambda(k) = run.lambda(1);
    resolutionResidual(k) = run.diagnostics.leadingRelativeResidual;
    if k==numel(elementsPerPeriod)
        reference = run;
    end
end
resolutionTable = table(elementsPerPeriod,resolutionLambda, ...
    resolutionResidual,'VariableNames', ...
    {'elementsPerPeriod','lambda1','relativeResidual'});

domainConvergence = abs(domainLambda(end)-domainLambda(end-1)) ...
    /max(1,abs(domainLambda(end)));
resolutionConvergence = abs( ...
    resolutionLambda(end)-resolutionLambda(end-1)) ...
    /max(1,abs(resolutionLambda(end)));
d = reference.diagnostics;

% A stable ideal-ballooning operator has a continuum edge at lambda=0.
% Its largest finite-Dirichlet eigenvalue then approaches zero from below
% as -constant/nPeriods^2 and is not a localized decaying mode. A positive
% discrete eigenvalue instead must converge and decay at the boundary.
inverseLengthSquared = 1./domainPeriods.^2;
fitMatrix = [ones(size(domainPeriods)),inverseLengthSquared];
fitCoefficient = fitMatrix\domainLambda;
lambdaInfinite = fitCoefficient(1);
fitPrediction = fitMatrix*fitCoefficient;
fitRelativeRms = sqrt(mean((fitPrediction-domainLambda).^2)) ...
    /max(max(abs(domainLambda)),eps);
continuumScaled = -domainLambda.*domainPeriods.^2;
continuumScalingChange = abs( ...
    continuumScaled(end)-continuumScaled(end-1)) ...
    /max(abs(continuumScaled(end)),eps);

if domainLambda(end)>0
    classification = 'unstable discrete mode';
    domainBehaviorPassed = domainConvergence<2e-3 ...
        && d.leadingEdgeAmplitudeRatio<1e-2;
else
    classification = 'stable continuum edge (lambda approaches 0-)';
    domainBehaviorPassed = all(domainLambda<0) ...
        && all(diff(domainLambda)>0) ...
        && continuumScalingChange<0.05 ...
        && abs(lambdaInfinite)<5e-4;
end

fprintf('Sparse/dense leading-eigenvalue error: %.3e\n',denseSparseError)
disp(domainTable)
disp(resolutionTable)
fprintf(['Reference mode: lambda1=%.9g, residual=%.3e, ' ...
         'edge amplitude=%.3e, parity error=%.3e\n'], ...
    reference.lambda(1),d.leadingRelativeResidual, ...
    d.leadingEdgeAmplitudeRatio,d.leadingParityError)
fprintf(['Infinite-domain interpretation: %s, lambda_inf=%.3e, ' ...
         'fit RMS=%.3e, N^2 scaling change=%.3e\n'], ...
    classification,lambdaInfinite,fitRelativeRms, ...
    continuumScalingChange)

assert(denseSparseError<1e-10, ...
    'Sparse and dense generalized eigenvalues disagree.');
assert(domainBehaviorPassed, ...
    'The field-line-domain behavior matches neither a discrete mode nor a stable continuum edge.');
assert(resolutionConvergence<2e-3, ...
    'The leading eigenvalue is not converged in angular resolution.');
assert(d.eigsFlag==0 && d.maximumRelativeResidual<1e-8, ...
    'The sparse generalized eigensolver did not converge accurately.');
assert(d.massCholeskyFlag==0, ...
    'The eigenvalue weight matrix is not positive definite.');
assert(d.boundaryAmplitude<1e-14, ...
    'The Dirichlet boundary condition was not imposed exactly.');
assert(d.rawOperatorSymmetryError<1e-13 && ...
       d.rawMassSymmetryError<1e-13, ...
    'The assembled weak-form matrices are not symmetric.');

report = struct();
report.localEquilibrium = localEq;
report.solution = reference;
report.domainConvergence = domainTable;
report.resolutionConvergence = resolutionTable;
report.denseSparseRelativeError = denseSparseError;
report.domainRelativeChange = domainConvergence;
report.resolutionRelativeChange = resolutionConvergence;
report.infiniteDomainLambdaEstimate = lambdaInfinite;
report.domainFitRelativeRms = fitRelativeRms;
report.continuumScalingChange = continuumScalingChange;
report.classification = classification;

fprintf('local ballooning eigensolver verification passed.\n')

if doPlot
    plot_local_ballooning_mode(reference);
end
end


function add_sturm_liouville_functions()
gsRoot = fileparts(fileparts(mfilename('fullpath')));
folder = fullfile(fileparts(gsRoot),'sturm_liouville_solver');
if ~isfolder(folder)
    error('GS:ballooning:MissingSturmLiouvilleFolder', ...
        'The sibling sturm_liouville_solver folder was not found.');
end
addpath(folder);
end
