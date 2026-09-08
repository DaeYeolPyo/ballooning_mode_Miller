function [psiNormalized, info] = normalize_flux( ...
    psi, psiAxis, psiBoundary, opts)
%NORMALIZE_FLUX Map magnetic flux to 0 at the axis and 1 at the LCFS.
%
%   PSIN = NORMALIZE_FLUX(PSI, PSIAXIS, PSIBOUNDARY) computes
%
%       PSIN = (PSI - PSIAXIS)/(PSIBOUNDARY - PSIAXIS).
%
%   The formula supports either flux orientation: PSIAXIS may be smaller or
%   larger than PSIBOUNDARY.
%
%   PSIN = NORMALIZE_FLUX(..., OPTS) accepts:
%       rangePolicy : 'extrapolate' (default), 'clip', or 'error'
%       tolerance   : normalized-coordinate tolerance (default 1e-10)
%
%   Under the 'error' policy, values outside [0,1] by more than tolerance
%   are rejected; roundoff-sized overshoots are clamped. INFO reports the
%   raw range and the number of values outside the physical interval.

narginchk(3, 4);

validateattributes(psi, {'numeric'}, {'real', 'finite'}, mfilename, 'psi');
validateattributes(psiAxis, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'psiAxis');
validateattributes(psiBoundary, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'psiBoundary');

span = psiBoundary-psiAxis;
fluxScale = max([1, abs(psiAxis), abs(psiBoundary)]);
if abs(span) <= 100*eps(fluxScale)
    error('GS:physics:DegenerateFluxRange', ...
        'psiAxis and psiBoundary must be numerically distinct.');
end

if nargin < 4 || isempty(opts)
    opts = struct();
end
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:physics:InvalidNormalizationOptions', ...
        'opts must be a scalar struct.');
end

defaults = struct('rangePolicy', 'extrapolate', 'tolerance', 1e-10);
unknown = setdiff(fieldnames(opts), fieldnames(defaults));
if ~isempty(unknown)
    error('GS:physics:UnknownNormalizationOption', ...
        'Unknown normalization option: %s', strjoin(unknown, ', '));
end

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts, name)
        opts.(name) = defaults.(name);
    end
end

validateattributes(opts.tolerance, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, ...
    mfilename, 'opts.tolerance');
policy = validate_policy(opts.rangePolicy);

raw = (psi-psiAxis)/span;
below = raw < 0;
above = raw > 1;
farBelow = raw < -opts.tolerance;
farAbove = raw > 1+opts.tolerance;

switch policy
    case 'extrapolate'
        psiNormalized = raw;
    case 'clip'
        psiNormalized = min(max(raw, 0), 1);
    case 'error'
        if any(farBelow(:)) || any(farAbove(:))
            error('GS:physics:FluxOutsidePlasma', ...
                ['Normalized flux lies outside [0,1]: raw range is ' ...
                 '[%.6g, %.6g].'], min(raw(:)), max(raw(:)));
        end
        psiNormalized = min(max(raw, 0), 1);
end

if nargout > 1
    info = struct();
    info.psiAxis = psiAxis;
    info.psiBoundary = psiBoundary;
    info.span = span;
    info.rangePolicy = policy;
    if isempty(raw)
        info.rawRange = [NaN, NaN];
    else
        info.rawRange = [min(raw(:)), max(raw(:))];
    end
    info.nBelow = nnz(below);
    info.nAbove = nnz(above);
    info.nOutsideTolerance = nnz(farBelow) + nnz(farAbove);
    info.nClipped = nnz(psiNormalized ~= raw);
end

end


function policy = validate_policy(value)
if isstring(value)
    if ~isscalar(value)
        error('GS:physics:InvalidRangePolicy', ...
            'rangePolicy must be a character vector or string scalar.');
    end
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:physics:InvalidRangePolicy', ...
        'rangePolicy must be a character vector or string scalar.');
end

policy = lower(strtrim(value));
allowed = {'extrapolate', 'clip', 'error'};
if ~ismember(policy, allowed)
    error('GS:physics:InvalidRangePolicy', ...
        'rangePolicy must be extrapolate, clip, or error.');
end
end
