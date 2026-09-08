function solution = solve_local_ballooning_mode(localEq,targetPsiN,opts)
%SOLVE_LOCAL_BALLOONING_MODE Solve the local infinite-field-line problem.
%
%   SOLUTION = SOLVE_LOCAL_BALLOONING_MODE(LOCALEQ,TARGETPSIN) builds the
%   regularized PEST coefficients, truncates the nominal infinite domain to
%   [-nPeriods*pi,nPeriods*pi], imposes X=0 at both ends, and solves
%
%       d/dtheta(g dX/dtheta) + c X = lambda f X.
%
%   Under this sign convention, lambda>0 is ballooning unstable. The P2
%   finite-element matrices are assembled by the existing sibling
%   sturm_liouville_solver/construct_matrix.m. Only the requested leading
%   eigenpairs are computed with a sparse symmetric generalized solve.
%
%   OPTS fields:
%       theta0              ballooning angle (default 0)
%       nPeriods            total number of 2*pi periods (default 6)
%       elementsPerPeriod   P2 elements per 2*pi (default 128)
%       numEigenvalues      leading eigenpairs retained (default 4)
%       poloidalModeCutoff  coefficient regularization cutoff (default 24)
%       eigenTolerance      eigs relative tolerance (default 1e-10)
%       maxIterations       eigs maximum iterations (default 2000)

narginchk(2,3);
if nargin<3 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts);
add_sturm_liouville_functions();
if exist('construct_matrix','file')~=2
    error('GS:ballooning:MissingMatrixAssembler', ...
        'The sibling construct_matrix.m function is unavailable.');
end

localBal = build_local_ballooning_coefficients( ...
    localEq,targetPsiN,struct( ...
        'theta0',opts.theta0, ...
        'nPeriods',opts.nPeriods, ...
        'poloidalModeCutoff',opts.poloidalModeCutoff));
coeff = localBal.coefficients;

nElements = opts.nPeriods*opts.elementsPerPeriod;
thetaGrid = linspace(-opts.nPeriods*pi,opts.nPeriods*pi, ...
    nElements+1).';
g = interp1(coeff.theta,coeff.g,thetaGrid,'pchip');
c = interp1(coeff.theta,coeff.c,thetaGrid,'pchip');
f = interp1(coeff.theta,coeff.f,thetaGrid,'pchip');
if any(~isfinite([g;c;f])) || any(g<=0) || any(f<=0)
    error('GS:ballooning:InvalidAssemblyCoefficient', ...
        'The resampled g and f must be finite and positive.');
end

[Gmat,Cmat,Fmat,thetaDof] = construct_matrix(thetaGrid,g,c,f);
rawA = Cmat-Gmat;
rawASymmetry = relative_symmetry_error(rawA);
rawFSymmetry = relative_symmetry_error(Fmat);
Amat = 0.5*(rawA+rawA.');
Fmat = 0.5*(Fmat+Fmat.');

ndof = numel(thetaDof);
free = (2:ndof-1).';
Afree = Amat(free,free);
Ffree = Fmat(free,free);
[~,massCholeskyFlag] = chol(Ffree);
if massCholeskyFlag~=0
    error('GS:ballooning:NonPositiveMassMatrix', ...
        'The Dirichlet mass matrix is not positive definite.');
end
if opts.numEigenvalues>=numel(free)
    error('GS:ballooning:TooManyEigenvalues', ...
        'numEigenvalues must be smaller than the free matrix dimension.');
end

% Rayleigh's principle gives lambda_max <= max(c/f). Use a shift safely
% above this upper bound, so sparse shift-invert returns the algebraically
% largest eigenvalues without competing with the unbounded negative
% field-line-bending spectrum.
cOverFUpperBound = max(c./f);
eigenShift = cOverFUpperBound+max(1,abs(cOverFUpperBound));
nFree = numel(free);
initialVector = ones(nFree,1)/sqrt(nFree);
eigsOptions = struct('tol',opts.eigenTolerance, ...
    'maxit',opts.maxIterations, ...
    'disp',0,'v0',initialVector);
[V,D,eigsFlag] = eigs(Afree,Ffree,opts.numEigenvalues, ...
    eigenShift,eigsOptions);
lambda = real(diag(D));
[lambda,order] = sort(lambda,'descend');
V = real(V(:,order));

X = zeros(ndof,opts.numEigenvalues);
residual = zeros(opts.numEigenvalues,1);
rayleigh = zeros(opts.numEigenvalues,1);
for modeIndex = 1:opts.numEigenvalues
    v = V(:,modeIndex);
    massNorm = sqrt(v.'*Ffree*v);
    v = v/massNorm;
    [~,peak] = max(abs(v));
    if v(peak)<0
        v = -v;
    end
    X(free,modeIndex) = v;
    Av = Afree*v;
    Fv = Ffree*v;
    residual(modeIndex) = norm(Av-lambda(modeIndex)*Fv) ...
        /max(norm(Av)+abs(lambda(modeIndex))*norm(Fv),eps);
    rayleigh(modeIndex) = (v.'*Av)/(v.'*Fv);
end

leadingMode = X(:,1);
thetaMax = max(abs(thetaDof));
edge = abs(thetaDof)>=0.9*thetaMax;
lumpedWeight = full(sum(Fmat,2));
totalEnergy = sum(lumpedWeight.*leadingMode.^2);
edgeEnergy = sum(lumpedWeight(edge).*leadingMode(edge).^2);
reflection = flipud(leadingMode);
parityError = min(norm(leadingMode-reflection), ...
    norm(leadingMode+reflection))/max(norm(leadingMode),eps);

diagnostics = struct();
diagnostics.eigsFlag = eigsFlag;
diagnostics.maximumRelativeResidual = max(residual);
diagnostics.leadingRelativeResidual = residual(1);
diagnostics.maximumRayleighError = max(abs(rayleigh-lambda)) ...
    /max(max(abs(lambda)),1);
diagnostics.rawOperatorSymmetryError = rawASymmetry;
diagnostics.rawMassSymmetryError = rawFSymmetry;
diagnostics.massCholeskyFlag = massCholeskyFlag;
diagnostics.rayleighUpperBound = cOverFUpperBound;
diagnostics.eigenShift = eigenShift;
diagnostics.boundaryAmplitude = max(abs( ...
    [X(1,:);X(end,:)]),[],'all');
diagnostics.leadingEdgeAmplitudeRatio = ...
    max(abs(leadingMode(edge)))/max(abs(leadingMode));
diagnostics.leadingEdgeEnergyFraction = edgeEnergy/totalEnergy;
diagnostics.leadingParityError = parityError;
diagnostics.leadingPeakTheta = thetaDof( ...
    find(abs(leadingMode)==max(abs(leadingMode)),1));

solution = struct();
solution.model = 'local PEST ideal ballooning eigenproblem';
solution.signConvention = 'lambda>0 is unstable';
if lambda(1)>0
    solution.finiteDomainClassification = ...
        'positive discrete candidate (unstable)';
else
    solution.finiteDomainClassification = ...
        'no positive finite-domain eigenvalue';
end
solution.targetPsiN = targetPsiN;
solution.thetaGrid = thetaGrid;
solution.thetaDof = thetaDof;
solution.g = g;
solution.c = c;
solution.f = f;
solution.lambda = lambda;
solution.X = X;
solution.leadingMode = leadingMode;
solution.localBallooning = localBal;
solution.options = opts;
solution.matrices = struct('G',Gmat,'C',Cmat,'F',Fmat, ...
    'A',Amat,'free',free);
solution.diagnostics = diagnostics;
end


function value = relative_symmetry_error(matrix)
value = norm(matrix-matrix.','fro')/max(norm(matrix,'fro'),eps);
end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:ballooning:InvalidEigenOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('theta0',0,'nPeriods',6, ...
    'elementsPerPeriod',128,'numEigenvalues',4, ...
    'poloidalModeCutoff',24,'eigenTolerance',1e-10, ...
    'maxIterations',2000);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:ballooning:UnknownEigenOption', ...
        'Unknown option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.theta0,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.theta0');
integerFields = {'nPeriods','elementsPerPeriod','numEigenvalues', ...
    'poloidalModeCutoff','maxIterations'};
for k = 1:numel(integerFields)
    name = integerFields{k};
    validateattributes(opts.(name),{'numeric'}, ...
        {'real','finite','scalar','integer','positive'},mfilename,name);
end
if opts.elementsPerPeriod<2*opts.poloidalModeCutoff+2
    error('GS:ballooning:AngularGridTooCoarse', ...
        ['elementsPerPeriod must exceed twice poloidalModeCutoff so the ' ...
         'retained coefficient spectrum is resolved.']);
end
validateattributes(opts.eigenTolerance,{'numeric'}, ...
    {'real','finite','scalar','positive','<',1}, ...
    mfilename,'opts.eigenTolerance');
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
