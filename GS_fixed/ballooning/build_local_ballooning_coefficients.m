function localBal = build_local_ballooning_coefficients( ...
        localEq,targetPsiN,opts)
%BUILD_LOCAL_BALLOONING_COEFFICIENTS Construct local PEST coefficients.
%
%   LOCALBAL = BUILD_LOCAL_BALLOONING_COEFFICIENTS(LOCALEQ,TARGETPSIN)
%   converts the SI equilibrium produced by build_local_pest_equilibrium to
%   the dimensionless convention used by the sibling ballooning routines,
%   constructs the PEST metric, and evaluates
%
%       d/dtheta (g dX/dtheta) + c X = lambda f X.
%
%   This stage constructs and verifies g, c, and f only. It does not solve
%   the eigenvalue problem.
%
%   Normalization:
%       Rbar   = R/aN,                  aN = magnetic-axis R
%       psibar = psi/(BN*aN^2),         BN = |F_axis|/aN
%       Fbar   = F/(BN*aN)
%       pbar   = mu0*p/BN^2.
%
%   OPTS fields:
%       theta0      ballooning angle theta_0 (default 0)
%       nPeriods    number of 2*pi periods in the returned domain
%                   (default 1, centered about theta=0)
%       poloidalModeCutoff
%                   highest retained Fourier mode in R, Z, lambda, and
%                   B^2 (default 24). This removes P1 mesh-scale ripple,
%                   not resolved equilibrium shaping.

narginchk(2,3);
if nargin<3 || isempty(opts)
    opts = struct();
end
validate_local_equilibrium(localEq);
opts = validate_options(opts);
add_external_functions();

required = {'compute_metrics','calculate_ballooning_coeffs'};
for k = 1:numel(required)
    if exist(required{k},'file')~=2
        error('GS:ballooning:MissingExternalFunction', ...
            'Required external function %s is unavailable.',required{k});
    end
end

psiN = localEq.map.psiN(:);
[targetError,targetIndex] = min(abs(psiN-targetPsiN));
tolerance = 100*eps(max(1,max(abs(psiN))));
if targetError>tolerance
    error('GS:ballooning:TargetSurfaceMissing', ...
        ['targetPsiN=%.16g is not present. Include it in the local ' ...
         'surface bundle; nearest available value is %.16g.'], ...
        targetPsiN,psiN(targetIndex));
end
if targetIndex<3 || targetIndex>numel(psiN)-2
    error('GS:ballooning:InsufficientTargetStencil', ...
        'The target surface needs two radial neighbors on each side.');
end

[mapRaw,eqfuncN,normalization] = normalize_for_ballooning(localEq);
if opts.poloidalModeCutoff>=floor(numel(mapRaw.theta)/2)
    error('GS:ballooning:ModeCutoffTooHigh', ...
        'poloidalModeCutoff must be below the angular Nyquist mode.');
end
mapN = regularize_poloidal_map(mapRaw,opts.poloidalModeCutoff);
metrics = compute_metrics(mapN,eqfuncN);
metrics.B2_unfiltered = metrics.B2;
metrics.B2 = truncate_fourier_rows( ...
    metrics.B2,opts.poloidalModeCutoff);
if any(metrics.B2(:)<=0)
    error('GS:ballooning:NonPositiveFilteredB2', ...
        'Fourier regularization produced B^2<=0. Reduce smoothing.');
end

nTheta = numel(mapN.theta);
nIntervals = opts.nPeriods*nTheta;
thetaExtent = linspace(-opts.nPeriods*pi,opts.nPeriods*pi, ...
    nIntervals+1);
coeff = calculate_ballooning_coeffs( ...
    mapN,metrics,targetPsiN, ...
    Theta0=opts.theta0,ThetaExt=thetaExtent);

diagnostics = calculate_diagnostics( ...
    localEq,mapRaw,mapN,eqfuncN,metrics,coeff,targetIndex,normalization);

physical = struct();
physical.psi = localEq.map.psi(targetIndex);
physical.q = localEq.profiles.q(targetIndex);
physical.qprime = localEq.profiles.qprime(targetIndex);
physical.pprime = localEq.map.pprime(targetIndex);
physical.sHat = localEq.profiles.sHat(targetIndex);
physical.alpha = localEq.profiles.alpha(targetIndex);

localBal = struct();
localBal.model = 'FEM-GS local PEST ballooning coefficients';
localBal.targetPsiN = psiN(targetIndex);
localBal.targetIndex = targetIndex;
localBal.normalization = normalization;
localBal.rawNormalizedMap = mapRaw;
localBal.normalizedMap = mapN;
localBal.normalizedEqfunc = eqfuncN;
localBal.metrics = metrics;
localBal.coefficients = coeff;
localBal.physicalReference = physical;
localBal.diagnostics = diagnostics;
localBal.options = opts;
localBal.externalFunctionsCalled = required;
end


function map = regularize_poloidal_map(map,modeCutoff)
map.R = truncate_fourier_rows(map.R,modeCutoff);
map.Z = truncate_fourier_rows(map.Z,modeCutoff);
map.lambda = truncate_fourier_rows(map.lambda,modeCutoff);
if isfield(map,'Jacobian')
    map.Jacobian = repmat(map.q_profile(:),1,numel(map.theta)) ...
        .*map.R.^2./repmat(map.F(:),1,numel(map.theta));
end
end


function filtered = truncate_fourier_rows(values,modeCutoff)
n = size(values,2);
if mod(n,2)==0
    modes = [0:n/2,-n/2+1:-1];
else
    modes = [0:(n-1)/2,-(n-1)/2:-1];
end
spectrum = fft(values,[],2);
spectrum(:,abs(modes)>modeCutoff) = 0;
filtered = real(ifft(spectrum,[],2));
end


function [mapN,eqfuncN,norms] = normalize_for_ballooning(localEq)
mu0 = 4*pi*1e-7;
aN = localEq.bundle.axis.R;
Faxis = localEq.eqfunc.F(0);
BN = abs(Faxis/aN);
if ~(isfinite(aN) && aN>0 && isfinite(BN) && BN>0)
    error('GS:ballooning:InvalidNormalization', ...
        'The axis radius and reference magnetic field must be positive.');
end

psiScale = BN*aN^2;
FScale = BN*aN;
pScale = BN^2/mu0;
pprimeScale = pScale/psiScale;

mapN = localEq.map;
mapN.R = localEq.map.R/aN;
mapN.Z = localEq.map.Z/aN;
mapN.psi = localEq.map.psi/psiScale;
mapN.F = localEq.map.F/FScale;
mapN.pprime = localEq.map.pprime/pprimeScale;

if isfield(mapN,'V')
    mapN.V = localEq.map.V/aN^3;
end
if isfield(mapN,'Vprime')
    mapN.Vprime = localEq.map.Vprime*BN/aN;
end
if isfield(mapN,'qprime')
    mapN.qprime = localEq.map.qprime*psiScale;
end
if isfield(mapN,'Jacobian')
    mapN.Jacobian = localEq.map.Jacobian*BN/aN;
end
if isfield(mapN,'psiEvaluated')
    mapN.psiEvaluated = localEq.map.psiEvaluated/psiScale;
end
if isfield(mapN,'psiR')
    mapN.psiR = localEq.map.psiR/(BN*aN);
end
if isfield(mapN,'psiZ')
    mapN.psiZ = localEq.map.psiZ/(BN*aN);
end
if isfield(mapN,'gradPsi')
    mapN.gradPsi = localEq.map.gradPsi/(BN*aN);
end
if isfield(mapN,'B2')
    mapN.B2 = localEq.map.B2/BN^2;
end
if isfield(mapN,'Bp')
    mapN.Bp = localEq.map.Bp/BN;
end
if isfield(mapN,'Bphi')
    mapN.Bphi = localEq.map.Bphi/BN;
end

eqfuncN = struct();
eqfuncN.psiR = @(Rbar,Zbar) localEq.eqfunc.psiR( ...
    aN.*Rbar,aN.*Zbar)/(BN*aN);
eqfuncN.psiZ = @(Rbar,Zbar) localEq.eqfunc.psiZ( ...
    aN.*Rbar,aN.*Zbar)/(BN*aN);

norms = struct('aN',aN,'BN',BN,'mu0',mu0, ...
    'psiScale',psiScale,'FScale',FScale,'pressureScale',pScale, ...
    'pprimeScale',pprimeScale, ...
    'definitions',struct( ...
        'length','Rbar=R/aN', ...
        'flux','psibar=psi/(BN*aN^2)', ...
        'field','Bbar=B/BN', ...
        'pressure','pbar=mu0*p/BN^2'));
end


function diagnostics = calculate_diagnostics( ...
        localEq,mapRaw,mapN,eqfuncN,metrics,coeff,k,norms)
Jgeometry = metrics.Jacobian(k,:);
Janalytic = mapN.q_profile(k).*mapN.R(k,:).^2/mapN.F(k);
gradPsi2 = eqfuncN.psiR(mapN.R(k,:),mapN.Z(k,:)).^2 ...
    +eqfuncN.psiZ(mapN.R(k,:),mapN.Z(k,:)).^2;

detgError = relative_rms(metrics.detg(k,:),Jgeometry.^2);
jacobianError = relative_rms(Jgeometry,Janalytic);
gradPsiError = relative_rms(metrics.g_contra.psipsi(k,:),gradPsi2);
B2MetricError = relative_rms( ...
    metrics.B2_metric(k,:),metrics.B2_direct(k,:));
filteredB2Change = relative_rms( ...
    metrics.B2(k,:),metrics.B2_unfiltered(k,:));
mapDisplacement = norms.aN*hypot( ...
    mapN.R-mapRaw.R,mapN.Z-mapRaw.Z);

expectedQprime = localEq.profiles.qprime(k)*norms.psiScale;
expectedPprime = localEq.map.pprime(k)/norms.pprimeScale;

% On the grid-aligned extended domain, f/g=J^2 is an exact algebraic
% consequence of the coefficient definitions. Evaluate the periodic J on
% the same theta grid to catch extension or indexing mistakes.
Jextended = periodic_sample(mapN.theta,Jgeometry,coeff.theta);
fgIdentityError = max(abs(coeff.f./coeff.g-Jextended.^2)) ...
    /max(max(abs(Jextended.^2)),eps);

diagnostics = struct();
diagnostics.allMetricJacobiansPositive = all(metrics.Jacobian(:)>0);
diagnostics.minimumMetricJacobian = min(metrics.Jacobian,[],'all');
diagnostics.maximumMetricJacobian = max(metrics.Jacobian,[],'all');
diagnostics.targetJacobianAnalyticRelativeRms = jacobianError;
diagnostics.targetMetricDeterminantRelativeRms = detgError;
diagnostics.targetGradPsiMetricRelativeRms = gradPsiError;
diagnostics.targetB2MetricRelativeRms = B2MetricError;
diagnostics.targetFilteredB2ChangeRelativeRms = filteredB2Change;
diagnostics.maximumMapRegularizationDisplacement = ...
    max(mapDisplacement,[],'all');
diagnostics.rmsMapRegularizationDisplacement = ...
    sqrt(mean(mapDisplacement.^2,'all'));
diagnostics.extendedFOverGIdentityRelativeError = fgIdentityError;
diagnostics.qprimeNormalizationRelativeError = ...
    abs(coeff.qprime-expectedQprime)/max(abs(expectedQprime),eps);
diagnostics.pprimeNormalizationRelativeError = ...
    abs(coeff.pprime-expectedPprime)/max(abs(expectedPprime),eps);
diagnostics.minimumK = min(coeff.K);
diagnostics.minimumG = min(coeff.g);
diagnostics.minimumF = min(coeff.f);
diagnostics.maximumAbsC = max(abs(coeff.c));
diagnostics.allCoefficientsFinite = all(isfinite( ...
    [coeff.g;coeff.c;coeff.f;coeff.K]));
end


function yq = periodic_sample(theta,y,thetaq)
theta = theta(:).';
y = y(:).';
theta0 = theta(1);
wrapped = mod(thetaq(:).'-theta0,2*pi)+theta0;
yq = interp1([theta,theta0+2*pi],[y,y(1)],wrapped,'pchip').';
end


function value = relative_rms(a,b)
value = sqrt(mean((a(:)-b(:)).^2))/max(sqrt(mean(b(:).^2)),eps);
end


function validate_local_equilibrium(localEq)
required = {'map','profiles','bundle','eqfunc'};
if ~isstruct(localEq) || ~isscalar(localEq) || ...
        ~all(isfield(localEq,required)) || ...
        ~all(isfield(localEq.eqfunc,{'F','psiR','psiZ'}))
    error('GS:ballooning:InvalidLocalEquilibrium', ...
        'localEq must be returned by build_local_pest_equilibrium.');
end
if numel(localEq.map.psiN)<5
    error('GS:ballooning:InsufficientRadialSurfaces', ...
        'At least five local flux surfaces are required.');
end
end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:ballooning:InvalidCoefficientOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('theta0',0,'nPeriods',1,'poloidalModeCutoff',24);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:ballooning:UnknownCoefficientOption', ...
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
validateattributes(opts.nPeriods,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'}, ...
    mfilename,'opts.nPeriods');
validateattributes(opts.poloidalModeCutoff,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',3}, ...
    mfilename,'opts.poloidalModeCutoff');
end


function add_external_functions()
gsRoot = fileparts(fileparts(mfilename('fullpath')));
parent = fileparts(gsRoot);
folders = {fullfile(parent,'equilibrium'), ...
           fullfile(parent,'ballooning_equation')};
for k = 1:numel(folders)
    if ~isfolder(folders{k})
        error('GS:ballooning:MissingExternalFolder', ...
            'Required sibling folder was not found: %s',folders{k});
    end
end
addpath(folders{:});
end
