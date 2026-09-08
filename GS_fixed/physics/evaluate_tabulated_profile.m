function values = evaluate_tabulated_profile(psiNormalized,spec,profileName)
%EVALUATE_TABULATED_PROFILE Interpolate a profile supplied on psi_N.

narginchk(2,3);
if nargin<3 || isempty(profileName)
    profileName = 'profile';
end
validateattributes(psiNormalized,{'numeric'},{'real','finite'}, ...
    mfilename,'psiNormalized');
if ~isstruct(spec) || ~isscalar(spec) || ...
        ~all(isfield(spec,{'psiN','values','method'}))
    error('GS:physics:InvalidTabulatedProfile', ...
        '%s must provide psiN, values, and method.',profileName);
end

knots = spec.psiN(:);
data = spec.values(:);
validateattributes(knots,{'numeric'}, ...
    {'real','finite','vector','nonempty'},mfilename,'spec.psiN');
validateattributes(data,{'numeric'}, ...
    {'real','finite','vector','numel',numel(knots)}, ...
    mfilename,'spec.values');
if numel(knots)<2 || any(diff(knots)<=0)
    error('GS:physics:InvalidTabulatedGrid', ...
        '%s psiN knots must be strictly increasing.',profileName);
end
gridTolerance = 1e-12;
if abs(knots(1))>gridTolerance || abs(knots(end)-1)>gridTolerance
    error('GS:physics:IncompleteTabulatedGrid', ...
        '%s psiN knots must span exactly from 0 to 1.',profileName);
end

method = normalize_method(spec.method,profileName);
values = interp1(knots,data,psiNormalized,method);
if ~isreal(values) || any(~isfinite(values(:)))
    error('GS:physics:InvalidTabulatedValues', ...
        '%s interpolation produced non-finite values.',profileName);
end
end


function method = normalize_method(value,profileName)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:physics:InvalidTabulatedMethod', ...
        '%s interpolation method must be text.',profileName);
end
method = lower(strtrim(value));
if ~ismember(method,{'linear','pchip','makima'})
    error('GS:physics:InvalidTabulatedMethod', ...
        '%s interpolation method must be linear, pchip, or makima.', ...
        profileName);
end
end
