function model = prepare_local_salpha_model(localEq,targetPsiN,opts)
%PREPARE_LOCAL_SALPHA_MODEL Freeze geometry and define s-alpha mappings.
%
%   The PEST geometry, q, F, and baseline radial B^2 derivative are frozen
%   at TARGETPSIN. Independent scan inputs are mapped through the same
%   volume definitions used by compute_contour_integrals:
%
%       qprime = sHat*q*Vprime/(2*V),
%       pprime = -alpha*(2*pi^2/Vprime)*sqrt(2*pi^2*Raxis/V).
%
%   Here qprime and pprime are derivatives with respect to dimensionless
%   flux, and pprime is a derivative of mu0*p/BN^2. The reference sHat and
%   alpha therefore reproduce the original equilibrium coefficients.

narginchk(2,3);
if nargin<3 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts);

base = build_local_ballooning_coefficients( ...
    localEq,targetPsiN,struct('theta0',0, ...
        'nPeriods',opts.nPeriods, ...
        'poloidalModeCutoff',opts.poloidalModeCutoff));
k = base.targetIndex;
map = base.normalizedMap;

q = map.q_profile(k);
V = map.V(k);
Vprime = map.Vprime(k);
Raxis = localEq.bundle.axis.R/base.normalization.aN;
if ~all(isfinite([q,V,Vprime,Raxis])) || ...
        q<=0 || V<=0 || Vprime<=0 || Raxis<=0
    error('GS:ballooning:InvalidSAlphaMapping', ...
        'q, V, Vprime, and normalized Raxis must be positive.');
end

qprimePerSHat = q*Vprime/(2*V);
pprimePerAlpha = -(2*pi^2/Vprime)*sqrt(2*pi^2*Raxis/V);

reference = struct();
reference.sHat = localEq.profiles.sHat(k);
reference.alpha = localEq.profiles.alpha(k);
reference.qprime = base.coefficients.qprime;
reference.pprime = base.coefficients.pprime;
reference.qprimeFromMapping = reference.sHat*qprimePerSHat;
reference.pprimeFromMapping = reference.alpha*pprimePerAlpha;
reference.qprimeRelativeError = abs( ...
    reference.qprimeFromMapping-reference.qprime) ...
    /max(abs(reference.qprime),eps);
reference.pprimeRelativeError = abs( ...
    reference.pprimeFromMapping-reference.pprime) ...
    /max(abs(reference.pprime),eps);

model = struct();
model.model = 'frozen-geometry local PEST s-alpha model';
model.targetPsiN = targetPsiN;
model.targetIndex = k;
model.localEquilibrium = localEq;
model.base = base;
model.q = q;
model.V = V;
model.Vprime = Vprime;
model.Raxis = Raxis;
model.qprimePerSHat = qprimePerSHat;
model.pprimePerAlpha = pprimePerAlpha;
model.reference = reference;
model.options = opts;
model.assumption = ['PEST geometry and baseline B2_psi are frozen; ' ...
    'qprime and pprime are varied independently.'];
end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:ballooning:InvalidSAlphaModelOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('nPeriods',6,'poloidalModeCutoff',24);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:ballooning:UnknownSAlphaModelOption', ...
        'Unknown option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.nPeriods,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'},mfilename,'nPeriods');
validateattributes(opts.poloidalModeCutoff,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',3}, ...
    mfilename,'poloidalModeCutoff');
end
