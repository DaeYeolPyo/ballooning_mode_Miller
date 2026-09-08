function [psi0, info] = initial_guess(meshData, K, opts)
%INITIAL_GUESS Build a smooth fixed-boundary Picard initial iterate.
%
%   PSI0 = INITIAL_GUESS(MESHDATA, K) solves a unit-source elliptic problem
%   on the actual LCFS, then scales its interior extremum to -1 while keeping
%   psi=0 exactly on the boundary. This produces a mesh-aware nested-looking
%   starting field without prescribing interior flux-surface geometry.
%
%   PSI0 = INITIAL_GUESS(..., OPTS) accepts:
%       boundaryValue : scalar LCFS flux (default 0)
%       axisValue     : desired interior extremum (default -1)
%
%   K may be empty, in which case assemble_stiffness(meshData) is called.

narginchk(1, 3);

if nargin < 2 || isempty(K)
    require_function('assemble_stiffness');
    K = assemble_stiffness(meshData);
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts);

required = {'nodes', 'boundaryNodes', 'interiorNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData, required))
    error('GS:solver:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
nNodes = size(meshData.nodes, 1);
if ~isnumeric(K) || ~isreal(K) || any(size(K) ~= [nNodes, nNodes])
    error('GS:solver:InvalidStiffnessMatrix', ...
        'K must be a real N-by-N stiffness matrix.');
end
if isempty(meshData.interiorNodes)
    error('GS:solver:NoInteriorNodes', ...
        'The mesh must contain at least one interior node.');
end
if abs(opts.axisValue-opts.boundaryValue) <= ...
        100*eps(max([1, abs(opts.axisValue), abs(opts.boundaryValue)]))
    error('GS:solver:DegenerateInitialFlux', ...
        'axisValue and boundaryValue must be numerically distinct.');
end

require_function('assemble_rhs');
require_function('apply_dirichlet');
unitLoad = assemble_rhs(meshData, 1.0);
[Kbc, rhsBc] = apply_dirichlet(meshData, K, unitLoad, 0.0);
bubble = Kbc\rhsBc;
bubble(meshData.boundaryNodes) = 0;

bubbleMaximum = max(bubble(meshData.interiorNodes));
if ~isfinite(bubbleMaximum) || bubbleMaximum <= 0
    error('GS:solver:InvalidBubbleFunction', ...
        'The unit-source initial-guess solve did not produce a positive interior.');
end

amplitude = opts.axisValue-opts.boundaryValue;
psi0 = opts.boundaryValue + amplitude*(bubble/bubbleMaximum);
psi0(meshData.boundaryNodes) = opts.boundaryValue;

if nargout > 1
    [axisValueActual, localIndex] = extremum_in_direction( ...
        psi0(meshData.interiorNodes), amplitude);
    info = struct();
    info.boundaryValue = opts.boundaryValue;
    info.requestedAxisValue = opts.axisValue;
    info.axisValue = axisValueActual;
    info.axisNode = meshData.interiorNodes(localIndex);
    info.axisPosition = meshData.nodes(info.axisNode,:);
    info.unitBubbleMaximum = bubbleMaximum;
end

end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:solver:InvalidInitialGuessOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('boundaryValue', 0.0, 'axisValue', -1.0);
unknown = setdiff(fieldnames(opts), fieldnames(defaults));
if ~isempty(unknown)
    error('GS:solver:UnknownInitialGuessOption', ...
        'Unknown initial-guess option: %s', strjoin(unknown, ', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts, names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.boundaryValue, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'opts.boundaryValue');
validateattributes(opts.axisValue, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'opts.axisValue');
end


function [value, index] = extremum_in_direction(values, amplitude)
if amplitude < 0
    [value, index] = min(values);
else
    [value, index] = max(values);
end
end


function require_function(name)
if exist(name, 'file') ~= 2
    error('GS:solver:MissingDependency', ...
        'Required function %s is not on the MATLAB path.', name);
end
end
