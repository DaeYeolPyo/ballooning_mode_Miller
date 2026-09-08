function [K, info] = assemble_stiffness(meshData)
%ASSEMBLE_STIFFNESS Assemble the P1 Grad-Shafranov stiffness matrix.
%
%   K = ASSEMBLE_STIFFNESS(MESHDATA) assembles
%
%       K_ij = integral_Omega (1/R) grad(N_i) dot grad(N_j) dR dZ
%
%   on the counterclockwise P1 triangles returned by generate_mesh.
%   The spatially varying 1/R coefficient is integrated with the symmetric
%   three-point, degree-two triangle rule.
%
%   [K, INFO] also returns assembly and symmetry diagnostics. Before a
%   Dirichlet condition is applied, K is positive semidefinite and the
%   constant vector is its expected null mode.

narginchk(1, 1);

[nodes, elements] = validate_mesh(meshData);
nNodes = size(nodes, 1);
nElements = size(elements, 1);

elementR = reshape(nodes(elements(:), 1), size(elements));
elementZ = reshape(nodes(elements(:), 2), size(elements));

R1 = elementR(:,1);
R2 = elementR(:,2);
R3 = elementR(:,3);
Z1 = elementZ(:,1);
Z2 = elementZ(:,2);
Z3 = elementZ(:,3);

twiceArea = (R2-R1).*(Z3-Z1) - (Z2-Z1).*(R3-R1);
if any(twiceArea <= 0)
    error('GS:fem:InvalidElementOrientation', ...
        ['All P1 elements must be nondegenerate and counterclockwise. ' ...
         'Use the mesh returned by generate_mesh.']);
end
elementArea = 0.5*twiceArea;

% Row e contains dN_i/dR or dN_i/dZ for i = 1,2,3 in element e.
gradR = [(Z2-Z3)./twiceArea, ...
         (Z3-Z1)./twiceArea, ...
         (Z1-Z2)./twiceArea];
gradZ = [(R3-R2)./twiceArea, ...
         (R1-R3)./twiceArea, ...
         (R2-R1)./twiceArea];

lambda = triangle_quadrature_points();
quadratureR = elementR*lambda.';
if any(quadratureR(:) <= 0)
    error('GS:fem:NonPositiveQuadratureRadius', ...
        'A triangle quadrature point reaches R <= 0.');
end

% Each quadrature point has reference weight 1/3 after multiplication by
% the physical triangle area.
integralInvR = (elementArea/3).*sum(1./quadratureR, 2);

rows = zeros(9*nElements, 1);
columns = zeros(9*nElements, 1);
values = zeros(9*nElements, 1);
entry = 0;
for i = 1:3
    for j = 1:3
        index = entry + (1:nElements);
        rows(index) = elements(:,i);
        columns(index) = elements(:,j);
        values(index) = integralInvR.*( ...
            gradR(:,i).*gradR(:,j) + gradZ(:,i).*gradZ(:,j));
        entry = entry + nElements;
    end
end

K = sparse(rows, columns, values, nNodes, nNodes);

if nargout > 1
    normK = norm(K, 'fro');
    info = struct();
    info.nNodes = nNodes;
    info.nElements = nElements;
    info.nNonzeros = nnz(K);
    info.quadrature = 'symmetric three-point';
    info.relativeSymmetryError = norm(K-K.', 'fro')/max(normK, eps);
    info.constantNullResidual = norm(K*ones(nNodes,1), inf) ...
                              / max(normK, eps);
end

end


function [nodes, elements] = validate_mesh(meshData)
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
if any(nodes(:,1) <= 0)
    error('GS:fem:NonPositiveNodeRadius', ...
        'All mesh nodes must lie at R > 0.');
end
end


function lambda = triangle_quadrature_points()
lambda = [2/3, 1/6, 1/6; ...
          1/6, 2/3, 1/6; ...
          1/6, 1/6, 2/3];
end
