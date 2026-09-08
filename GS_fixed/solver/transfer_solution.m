function [psiNew, info] = transfer_solution( ...
    oldMesh, oldPsi, newMesh, boundaryValue)
%TRANSFER_SOLUTION Transfer a P1 solution between continuation meshes.
%
%   PSINEW = TRANSFER_SOLUTION(OLDMESH, OLDPSI, NEWMESH, BOUNDARYVALUE)
%   uses barycentric interpolation inside the old constrained triangulation.
%   New-domain points outside the old plasma use nearest-node extrapolation.
%   The new LCFS is then set exactly to BOUNDARYVALUE.

narginchk(4, 4);

validate_mesh(oldMesh, 'oldMesh');
validate_mesh(newMesh, 'newMesh');
validateattributes(oldPsi, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', size(oldMesh.nodes,1)}, ...
    mfilename, 'oldPsi');
validateattributes(boundaryValue, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'boundaryValue');
oldPsi = oldPsi(:);

oldTriangulation = triangulation(oldMesh.elements, oldMesh.nodes);
triangleIndex = pointLocation(oldTriangulation, newMesh.nodes);
inside = ~isnan(triangleIndex);
psiNew = NaN(size(newMesh.nodes,1), 1);

if any(inside)
    barycentric = cartesianToBarycentric(oldTriangulation, ...
        triangleIndex(inside), newMesh.nodes(inside,:));
    oldElements = oldMesh.elements(triangleIndex(inside),:);
    oldElementValues = reshape(oldPsi(oldElements(:)), size(oldElements));
    psiNew(inside) = sum(barycentric.*oldElementValues, 2);
end

outside = ~inside;
if any(outside)
    nearestInterpolant = scatteredInterpolant( ...
        oldMesh.nodes(:,1), oldMesh.nodes(:,2), oldPsi, ...
        'nearest', 'nearest');
    psiNew(outside) = nearestInterpolant( ...
        newMesh.nodes(outside,1), newMesh.nodes(outside,2));
end

psiNew(newMesh.boundaryNodes) = boundaryValue;
if any(~isfinite(psiNew))
    error('GS:solver:SolutionTransferFailed', ...
        'Solution transfer produced a non-finite value.');
end

if nargout > 1
    info = struct();
    info.nInterpolated = nnz(inside);
    info.nExtrapolated = nnz(outside);
    info.extrapolatedFraction = nnz(outside)/numel(outside);
    info.boundaryValue = boundaryValue;
end

end


function validate_mesh(meshData, name)
required = {'nodes', 'elements', 'boundaryNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData, required))
    error('GS:solver:InvalidTransferMesh', ...
        '%s must be a mesh returned by generate_mesh.', name);
end
end
