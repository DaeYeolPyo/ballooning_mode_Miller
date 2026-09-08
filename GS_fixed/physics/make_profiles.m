function profiles = make_profiles(overrides)
%MAKE_PROFILES Create a validated Grad-Shafranov profile configuration.
%
%   PROFILES = MAKE_PROFILES() returns a zero-source configuration with SI
%   mu0 and the normalized-flux convention psiAxis=-1, psiBoundary=0.
%
%   PROFILES = MAKE_PROFILES(OVERRIDES) replaces top-level pprime and
%   FFprime specifications and merges individual fields of OVERRIDES.flux.
%   Supported top-level fields are mu0, flux, pprime, and FFprime.

narginchk(0, 1);

profiles = struct();
profiles.mu0 = 4*pi*1e-7;
profiles.flux = struct( ...
    'psiAxis', -1.0, ...
    'psiBoundary', 0.0, ...
    'rangePolicy', 'error', ...
    'tolerance', 1e-10);
profiles.pprime = struct('type', 'constant', 'value', 0.0);
profiles.FFprime = struct('type', 'constant', 'value', 0.0);

if nargin == 0 || isempty(overrides)
    overrides = struct();
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('GS:physics:InvalidProfileOverrides', ...
        'overrides must be a scalar struct.');
end

allowed = fieldnames(profiles);
unknown = setdiff(fieldnames(overrides), allowed);
if ~isempty(unknown)
    error('GS:physics:UnknownProfileOverride', ...
        'Unknown profile override: %s', strjoin(unknown, ', '));
end

if isfield(overrides, 'mu0')
    profiles.mu0 = overrides.mu0;
end
if isfield(overrides, 'pprime')
    profiles.pprime = overrides.pprime;
end
if isfield(overrides, 'FFprime')
    profiles.FFprime = overrides.FFprime;
end
if isfield(overrides, 'flux')
    if ~isstruct(overrides.flux) || ~isscalar(overrides.flux)
        error('GS:physics:InvalidFluxOverride', ...
            'overrides.flux must be a scalar struct.');
    end
    fluxFields = fieldnames(profiles.flux);
    unknownFlux = setdiff(fieldnames(overrides.flux), fluxFields);
    if ~isempty(unknownFlux)
        error('GS:physics:UnknownFluxOverride', ...
            'Unknown flux override: %s', strjoin(unknownFlux, ', '));
    end
    names = fieldnames(overrides.flux);
    for k = 1:numel(names)
        profiles.flux.(names{k}) = overrides.flux.(names{k});
    end
end

validateattributes(profiles.mu0, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'profiles.mu0');

probePsi = [profiles.flux.psiAxis; ...
            0.5*(profiles.flux.psiAxis+profiles.flux.psiBoundary); ...
            profiles.flux.psiBoundary];
normalize_flux(probePsi, profiles.flux.psiAxis, ...
    profiles.flux.psiBoundary, struct( ...
    'rangePolicy', profiles.flux.rangePolicy, ...
    'tolerance', profiles.flux.tolerance));
pprime(probePsi, profiles.pprime, profiles.flux);
FFprime(probePsi, profiles.FFprime, profiles.flux);

end
