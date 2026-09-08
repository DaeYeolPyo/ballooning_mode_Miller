function [rhs, info] = assemble_rhs(meshData, source, psi)
%ASSEMBLE_RHS Assemble a P1 load vector for the GS weak source.
%
%   RHS = ASSEMBLE_RHS(MESHDATA, SOURCE) assembles
%
%       rhs_i = integral_Omega q(R,Z) N_i dR dZ.
%
%   SOURCE may be:
%       * a scalar constant q,
%       * an N-node vector interpolated with the P1 basis, or
%       * a vectorized function handle q = source(R,Z).
%
%   RHS = ASSEMBLE_RHS(MESHDATA, SOURCE, PSI) also supports a vectorized
%   nonlinear callback q = source(R,Z,psi), where PSI is a nodal vector and
%   is interpolated to the same three quadrature points.
%
%   IMPORTANT SOURCE CONVENTION
%   ---------------------------
%   For the strong Grad-Shafranov equation
%
%       -Delta* psi = mu0*R^2*p'(psi) + FF'(psi),
%
%   pass the divided weak-form source
%
%       q = mu0*R*p'(psi) + FF'(psi)/R.
%
%   [RHS, INFO] returns quadrature and integral diagnostics.

narginchk(2, 3);

[nodes, elements, elementArea] = validate_mesh(meshData);
nNodes = size(nodes, 1);
nElements = size(elements, 1);

if nargin < 3
    psi = [];
elseif ~isempty(psi)
    validateattributes(psi, {'numeric'}, ...
        {'real', 'finite', 'vector', 'numel', nNodes}, mfilename, 'psi');
    psi = psi(:);
end

lambda = triangle_quadrature_points();
elementR = reshape(nodes(elements(:), 1), size(elements));
elementZ = reshape(nodes(elements(:), 2), size(elements));
quadratureR = elementR*lambda.';
quadratureZ = elementZ*lambda.';

if isnumeric(source)
    if isscalar(source)
        sourceAtQuadrature = repmat(source, nElements, 3);
    elseif isvector(source) && numel(source) == nNodes
        source = source(:);
        elementSource = reshape(source(elements(:)), size(elements));
        sourceAtQuadrature = elementSource*lambda.';
    else
        error('GS:fem:InvalidNumericSource', ...
            'Numeric source must be a scalar or one value per mesh node.');
    end
elseif isa(source, 'function_handle')
    if isempty(psi)
        quadraturePsi = [];
    else
        elementPsi = reshape(psi(elements(:)), size(elements));
        quadraturePsi = elementPsi*lambda.';
    end
    sourceAtQuadrature = evaluate_source( ...
        source, quadratureR, quadratureZ, quadraturePsi);
else
    error('GS:fem:InvalidSource', ...
        'source must be numeric or a function handle.');
end

if ~isnumeric(sourceAtQuadrature) || ~isreal(sourceAtQuadrature)
    error('GS:fem:InvalidSourceOutput', ...
        'The evaluated source must be a real numeric array.');
end
if isscalar(sourceAtQuadrature)
    sourceAtQuadrature = repmat(sourceAtQuadrature, nElements, 3);
elseif numel(sourceAtQuadrature) == 3*nElements
    sourceAtQuadrature = reshape(sourceAtQuadrature, nElements, 3);
else
    error('GS:fem:InvalidSourceOutputSize', ...
        ['The source callback must return a scalar or an array with the same ' ...
         'size as its R and Z inputs.']);
end
if any(~isfinite(sourceAtQuadrature(:)))
    error('GS:fem:NonFiniteSource', ...
        'The evaluated source contains NaN or Inf.');
end

% localLoad(e,i) = A_e/3 * sum_q q_e(q)*N_i(q).
localLoad = bsxfun(@times, elementArea/3, sourceAtQuadrature*lambda);
rhs = accumarray(elements(:), localLoad(:), [nNodes, 1], @sum, 0);

if nargout > 1
    info = struct();
    info.nNodes = nNodes;
    info.nElements = nElements;
    info.quadrature = 'symmetric three-point';
    info.integratedSource = sum(bsxfun(@times, ...
        elementArea/3, sourceAtQuadrature), 'all');
    info.loadSum = sum(rhs);
end

end


function values = evaluate_source(source, R, Z, psi)
nInputs = nargin(source);

if nInputs == 0
    values = source();
elseif nInputs == 1
    values = source([R(:), Z(:)]);
    if ~isscalar(values)
        values = reshape(values, size(R));
    end
elseif nInputs == 2
    values = source(R, Z);
elseif nInputs >= 3
    if isempty(psi)
        error('GS:fem:MissingPsi', ...
            ['This source callback accepts psi, but no nodal psi vector was ' ...
             'provided to assemble_rhs.']);
    end
    values = source(R, Z, psi);
else
    % A varargin callback has no fixed arity. Prefer the nonlinear form when
    % psi is available and otherwise use the two-coordinate form.
    if isempty(psi)
        values = source(R, Z);
    else
        values = source(R, Z, psi);
    end
end
end


function [nodes, elements, elementArea] = validate_mesh(meshData)
required = {'nodes', 'elements'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData, required))
    error('GS:fem:InvalidMesh', ...
        'meshData must be a scalar struct returned by generate_mesh.');
end

nodes = meshData.nodes;
elements = meshData.elements;
if ~isnumeric(nodes) || size(nodes,2) ~= 2 || isempty(nodes) || ...
        any(~isfinite(nodes(:)))
    error('GS:fem:InvalidNodes', ...
        'meshData.nodes must be a finite N-by-2 numeric array.');
end
if ~isnumeric(elements) || size(elements,2) ~= 3 || isempty(elements) || ...
        any(~isfinite(elements(:))) || any(elements(:) ~= round(elements(:))) || ...
        any(elements(:) < 1) || any(elements(:) > size(nodes,1))
    error('GS:fem:InvalidElements', ...
        'meshData.elements must contain valid one-based triangle indices.');
end

p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
twiceArea = (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) ...
          - (p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1));
if any(twiceArea <= 0)
    error('GS:fem:InvalidElementOrientation', ...
        'All mesh triangles must be counterclockwise and nondegenerate.');
end
elementArea = 0.5*twiceArea;
end


function lambda = triangle_quadrature_points()
lambda = [2/3, 1/6, 1/6; ...
          1/6, 2/3, 1/6; ...
          1/6, 1/6, 2/3];
end
