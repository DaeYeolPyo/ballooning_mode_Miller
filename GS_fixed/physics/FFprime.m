function values = FFprime(psi, spec, fluxSpec)
%FFPRIME Evaluate F(psi)*dF/dpsi for the GS source.
%
%   VALUES = FFPRIME(PSI, C) returns the constant C.
%
%   VALUES = FFPRIME(PSI, SPEC, FLUXSPEC) supports:
%
%   constant
%       SPEC.value
%
%   normalized-power
%       FF'(psi) = axisValue*(1-psi_N)^exponent
%       SPEC.axisValue, SPEC.exponent
%
%   tabulated
%       SPEC.psiN and SPEC.values contain FF' on [0,1]. SPEC.method is
%       linear, pchip, or makima. Optional SPEC.FValues stores F itself.
%
%   F-power
%       F(psi) = boundaryF
%              + (axisF-boundaryF)*(1-psi_N)^exponent
%       and this function returns F*dF/dpsi analytically.
%
%   FLUXSPEC contains psiAxis, psiBoundary, and optionally rangePolicy and
%   tolerance. It is unnecessary for a constant profile.

narginchk(2, 3);
validateattributes(psi, {'numeric'}, {'real', 'finite'}, mfilename, 'psi');
if nargin < 3
    fluxSpec = [];
end

if isnumeric(spec)
    validateattributes(spec, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'spec');
    values = repmat(spec, size(psi));
    return
end

if ~isstruct(spec) || ~isscalar(spec) || ~isfield(spec, 'type')
    error('GS:physics:InvalidFFPrimeSpec', ...
        'spec must be a numeric constant or a scalar struct with type.');
end
type = profile_type(spec.type);

switch type
    case 'constant'
        require_fields(spec, {'type', 'value'}, 'FFprime');
        validateattributes(spec.value, {'numeric'}, ...
            {'real', 'finite', 'scalar'}, mfilename, 'spec.value');
        values = repmat(spec.value, size(psi));

    case 'normalized-power'
        require_fields(spec, {'type', 'axisValue', 'exponent'}, 'FFprime');
        validateattributes(spec.axisValue, {'numeric'}, ...
            {'real', 'finite', 'scalar'}, mfilename, 'spec.axisValue');
        validateattributes(spec.exponent, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'nonnegative'}, ...
            mfilename, 'spec.exponent');
        [psiNormalized, ~] = normalized_coordinate(psi, fluxSpec);
        values = spec.axisValue*(1-psiNormalized).^spec.exponent;

    case 'f-power'
        allowed = {'type', 'axisF', 'boundaryF', 'exponent'};
        require_fields(spec, allowed, 'FFprime');
        validateattributes(spec.axisF, {'numeric'}, ...
            {'real', 'finite', 'scalar'}, mfilename, 'spec.axisF');
        validateattributes(spec.boundaryF, {'numeric'}, ...
            {'real', 'finite', 'scalar'}, mfilename, 'spec.boundaryF');
        validateattributes(spec.exponent, {'numeric'}, ...
            {'real', 'finite', 'scalar', '>=', 1}, ...
            mfilename, 'spec.exponent');
        [psiNormalized, span] = normalized_coordinate(psi, fluxSpec);
        deltaF = spec.axisF-spec.boundaryF;
        F = spec.boundaryF + deltaF*(1-psiNormalized).^spec.exponent;
        derivativeF = -(deltaF*spec.exponent/span) ...
                    *(1-psiNormalized).^(spec.exponent-1);
        values = F.*derivativeF;

    case 'tabulated'
        required = {'type','psiN','values','method'};
        allowed = [required,{'FValues'}];
        require_required_fields(spec,required,allowed,'FFprime');
        if isfield(spec,'FValues')
            validateattributes(spec.FValues,{'numeric'}, ...
                {'real','finite','vector','numel',numel(spec.psiN)}, ...
                mfilename,'spec.FValues');
        end
        [psiNormalized,~] = normalized_coordinate(psi,fluxSpec);
        values = evaluate_tabulated_profile( ...
            psiNormalized,spec,'FFprime');

    otherwise
        error('GS:physics:UnknownFFPrimeType', ...
            'Unknown FFprime profile type: %s', type);
end

if ~isreal(values) || any(~isfinite(values(:)))
    error('GS:physics:InvalidFFPrimeValues', ...
        ['FFprime produced complex or non-finite values. Check the flux range ' ...
         'policy and profile exponent.']);
end

end


function [psiNormalized, span] = normalized_coordinate(psi, fluxSpec)
if nargin < 2 || ~isstruct(fluxSpec) || ~isscalar(fluxSpec)
    error('GS:physics:MissingFluxSpec', ...
        'A scalar fluxSpec is required for a normalized FFprime profile.');
end
required = {'psiAxis', 'psiBoundary'};
allowed = [required, {'rangePolicy', 'tolerance'}];
unknown = setdiff(fieldnames(fluxSpec), allowed);
if ~all(isfield(fluxSpec, required)) || ~isempty(unknown)
    error('GS:physics:InvalidFluxSpec', ...
        ['fluxSpec must contain psiAxis and psiBoundary; optional fields are ' ...
         'rangePolicy and tolerance.']);
end

opts = struct();
if isfield(fluxSpec, 'rangePolicy')
    opts.rangePolicy = fluxSpec.rangePolicy;
else
    opts.rangePolicy = 'error';
end
if isfield(fluxSpec, 'tolerance')
    opts.tolerance = fluxSpec.tolerance;
end
psiNormalized = normalize_flux( ...
    psi, fluxSpec.psiAxis, fluxSpec.psiBoundary, opts);
span = fluxSpec.psiBoundary-fluxSpec.psiAxis;
end


function type = profile_type(value)
if isstring(value)
    if ~isscalar(value)
        error('GS:physics:InvalidProfileType', ...
            'Profile type must be a character vector or string scalar.');
    end
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:physics:InvalidProfileType', ...
        'Profile type must be a character vector or string scalar.');
end
type = lower(strrep(strtrim(value), '_', '-'));
end


function require_fields(spec, allowed, profileName)
unknown = setdiff(fieldnames(spec), allowed);
missing = setdiff(allowed, fieldnames(spec));
if ~isempty(unknown) || ~isempty(missing)
    error('GS:physics:InvalidProfileFields', ...
        'Invalid %s fields. Missing: [%s]. Unknown: [%s].', ...
        profileName, strjoin(missing, ', '), strjoin(unknown, ', '));
end
end


function require_required_fields(spec,required,allowed,profileName)
unknown = setdiff(fieldnames(spec),allowed);
missing = setdiff(required,fieldnames(spec));
if ~isempty(unknown) || ~isempty(missing)
    error('GS:physics:InvalidProfileFields', ...
        'Invalid %s fields. Missing: [%s]. Unknown: [%s].', ...
        profileName,strjoin(missing,', '),strjoin(unknown,', '));
end
end
