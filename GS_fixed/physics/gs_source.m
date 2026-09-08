function [source, components] = gs_source(R, Z, psi, profiles)
%GS_SOURCE Evaluate the divided weak-form Grad-Shafranov source.
%
%   Q = GS_SOURCE(R,Z,PSI,PROFILES) returns
%
%       Q = mu0*R*p'(PSI) + FF'(PSI)/R,
%
%   which is the source expected by fem/assemble_rhs.m. R, Z, and PSI must
%   be real arrays of the same size; scalar PSI is expanded. Every R must be
%   positive. PROFILES is created by make_profiles.
%
%   [Q, COMPONENTS] also returns the evaluated pprime, FFprime, pressure
%   contribution, and poloidal-current contribution.

narginchk(4, 4);

validateattributes(R, {'numeric'}, {'real', 'finite'}, mfilename, 'R');
validateattributes(Z, {'numeric'}, {'real', 'finite'}, mfilename, 'Z');
validateattributes(psi, {'numeric'}, {'real', 'finite'}, mfilename, 'psi');
if ~isequal(size(R), size(Z))
    error('GS:physics:CoordinateSizeMismatch', ...
        'R and Z must have the same size.');
end
if isscalar(psi)
    psi = repmat(psi, size(R));
elseif ~isequal(size(psi), size(R))
    error('GS:physics:FluxSizeMismatch', ...
        'psi must be scalar or have the same size as R and Z.');
end
if any(R(:) <= 0)
    error('GS:physics:NonPositiveRadius', ...
        'The Grad-Shafranov source requires R > 0.');
end

required = {'mu0', 'pprime', 'FFprime'};
if ~isstruct(profiles) || ~isscalar(profiles) || ...
        ~all(isfield(profiles, required))
    error('GS:physics:InvalidProfiles', ...
        'profiles must be a scalar configuration from make_profiles.');
end
validateattributes(profiles.mu0, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'profiles.mu0');

if isfield(profiles, 'flux')
    fluxSpec = profiles.flux;
else
    fluxSpec = [];
end

pprimeValue = pprime(psi, profiles.pprime, fluxSpec);
FFprimeValue = FFprime(psi, profiles.FFprime, fluxSpec);
pressureSource = profiles.mu0*R.*pprimeValue;
currentSource = FFprimeValue./R;
source = pressureSource+currentSource;

if any(~isfinite(source(:))) || ~isreal(source)
    error('GS:physics:InvalidGSSource', ...
        'The evaluated Grad-Shafranov weak source is non-finite or complex.');
end

if nargout > 1
    components = struct();
    components.pprime = pprimeValue;
    components.FFprime = FFprimeValue;
    components.pressure = pressureSource;
    components.poloidalCurrent = currentSource;
end

end
