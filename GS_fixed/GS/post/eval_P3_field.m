function [value, gradient, elementID] = eval_P3_field( ...
    P1, P3, nodalValue, Xq)
%EVAL_P3_FIELD Evaluate a scalar P3 FEM field and its physical gradient.
%
%   [value, gradient, elementID] = EVAL_P3_FIELD( ...
%       P1, P3, nodalValue, Xq)
%
% INPUT
%   P1
%       triangulation constructed from P3.P1elements/P3.P1points.
%   P3
%       P3 mesh structure.
%   nodalValue
%       One scalar value per global P3 degree of freedom.
%   Xq
%       Nq-by-2 Cartesian query points [R,Z].
%
% OUTPUT
%   value
%       Nq-by-1 interpolated field values.
%   gradient
%       Nq-by-2 physical derivatives [d/dR,d/dZ].
%   elementID
%       Containing P1 triangle for each query point.
%
% Points outside the triangulation return NaN in every output.  The P3
% solution is C0 conforming, so its elementwise gradient can be
% discontinuous on an element boundary.  pointLocation selects one of the
% adjacent elements for a query point lying exactly on such a boundary.

    narginchk(4, 4);

    requiredFields = {'points', 'elements', 'P1points', 'P1elements'};
    for k = 1:numel(requiredFields)
        name = requiredFields{k};
        if ~isfield(P3, name)
            error('GS:EvalP3Field:MissingMeshField', ...
                'P3.%s is required.', name);
        end
    end

    if ~isnumeric(Xq) || size(Xq, 2) ~= 2 || ...
            any(~isfinite(Xq), 'all')
        error('GS:EvalP3Field:InvalidQuery', ...
            'Xq must be a finite Nq-by-2 numeric array.');
    end

    nodalValue = nodalValue(:);
    nDof = size(P3.points, 1);

    if numel(nodalValue) ~= nDof || any(~isfinite(nodalValue))
        error('GS:EvalP3Field:InvalidNodalValue', ...
            ['nodalValue must contain one finite value per ', ...
             'global P3 degree of freedom.']);
    end

    nQuery = size(Xq, 1);

    value = nan(nQuery, 1);
    gradient = nan(nQuery, 2);

    elementID = pointLocation(P1, Xq);
    inside = find(~isnan(elementID));

    for k = inside(:).'
        e = elementID(k);

        lambda = cartesianToBarycentric(P1, e, Xq(k, :));
        xi = lambda(2);
        eta = lambda(3);

        vertexIDs = P3.P1elements(e, :);
        Ve = P3.P1points(vertexIDs, :);

        [N, dNphysical] = eval_P3_element(xi, eta, Ve);

        solutionIDs = P3.elements(e, :);
        localValue = nodalValue(solutionIDs);

        value(k) = N*localValue;
        gradient(k, :) = (dNphysical*localValue).';
    end
end
