function [Kbc, rhsBc, bc] = apply_dirichlet( ...
    meshData, K, rhs, boundaryValues)
%APPLY_DIRICHLET Impose fixed-LCFS values by symmetric elimination.
%
%   [KBC, RHSBC, BC] = APPLY_DIRICHLET(MESHDATA, K, RHS, G) imposes psi=G
%   on meshData.boundaryNodes. G may be:
%       * one scalar,
%       * one value per boundary node,
%       * one value per mesh node, or
%       * a vectorized function handle G(R,Z).
%
%   Boundary rows and columns are replaced by identity rows/columns after
%   their known contribution is moved to the right-hand side. Thus symmetry
%   and positive definiteness are preserved. BC contains boundaryNodes,
%   freeNodes, and prescribedValues for reconstruction and diagnostics.

narginchk(4, 4);

required = {'nodes', 'boundaryNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData, required))
    error('GS:fem:InvalidMesh', ...
        ['meshData must contain nodes and boundaryNodes from ' ...
         'generate_mesh.']);
end

nodes = meshData.nodes;
nNodes = size(nodes, 1);
if ~isnumeric(nodes) || size(nodes,2) ~= 2 || isempty(nodes) || ...
        any(~isfinite(nodes(:)))
    error('GS:fem:InvalidNodes', ...
        'meshData.nodes must be a finite N-by-2 numeric array.');
end
if ~isnumeric(K) || ~isreal(K) || ~ismatrix(K) || ...
        any(size(K) ~= [nNodes, nNodes]) || any(~isfinite(nonzeros(K)))
    error('GS:fem:InvalidMatrix', ...
        'K must be a finite N-by-N numeric matrix.');
end
if ~isnumeric(rhs) || ~isreal(rhs) || ~isvector(rhs) || ...
        numel(rhs) ~= nNodes || any(~isfinite(rhs(:)))
    error('GS:fem:InvalidRightHandSide', ...
        'rhs must be a finite real vector with one entry per mesh node.');
end
rhs = rhs(:);

boundaryNodes = meshData.boundaryNodes(:);
if isempty(boundaryNodes) || ...
        any(boundaryNodes ~= round(boundaryNodes)) || ...
        any(boundaryNodes < 1) || any(boundaryNodes > nNodes) || ...
        numel(unique(boundaryNodes)) ~= numel(boundaryNodes)
    error('GS:fem:InvalidBoundaryNodes', ...
        'meshData.boundaryNodes must contain unique valid node indices.');
end
nBoundary = numel(boundaryNodes);

if isnumeric(boundaryValues)
    if isscalar(boundaryValues)
        prescribedValues = repmat(boundaryValues, nBoundary, 1);
    elseif isvector(boundaryValues) && numel(boundaryValues) == nBoundary
        prescribedValues = boundaryValues(:);
    elseif isvector(boundaryValues) && numel(boundaryValues) == nNodes
        boundaryValues = boundaryValues(:);
        prescribedValues = boundaryValues(boundaryNodes);
    else
        error('GS:fem:InvalidBoundaryValues', ...
            ['Numeric boundaryValues must be scalar, one value per boundary ' ...
             'node, or one value per mesh node.']);
    end
elseif isa(boundaryValues, 'function_handle')
    R = nodes(boundaryNodes,1);
    Z = nodes(boundaryNodes,2);
    nInputs = nargin(boundaryValues);
    if nInputs == 1
        prescribedValues = boundaryValues([R, Z]);
    else
        prescribedValues = boundaryValues(R, Z);
    end
    if isscalar(prescribedValues)
        prescribedValues = repmat(prescribedValues, nBoundary, 1);
    else
        prescribedValues = prescribedValues(:);
    end
    if numel(prescribedValues) ~= nBoundary
        error('GS:fem:InvalidBoundaryFunctionOutput', ...
            ['The boundary callback must return a scalar or one value per ' ...
             'boundary node.']);
    end
else
    error('GS:fem:InvalidBoundaryValues', ...
        'boundaryValues must be numeric or a function handle.');
end

if ~isnumeric(prescribedValues) || ~isreal(prescribedValues) || ...
        any(~isfinite(prescribedValues))
    error('GS:fem:NonFiniteBoundaryValues', ...
        'Prescribed boundary values must be finite and real.');
end

rhsBc = rhs - K(:,boundaryNodes)*prescribedValues;
Kbc = K;
Kbc(:,boundaryNodes) = 0;
Kbc(boundaryNodes,:) = 0;
Kbc = Kbc + sparse(boundaryNodes, boundaryNodes, ...
    ones(nBoundary,1), nNodes, nNodes);
rhsBc(boundaryNodes) = prescribedValues;

isBoundary = false(nNodes, 1);
isBoundary(boundaryNodes) = true;
bc = struct();
bc.boundaryNodes = boundaryNodes;
bc.freeNodes = find(~isBoundary);
bc.prescribedValues = prescribedValues;

end
