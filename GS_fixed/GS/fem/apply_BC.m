function [KII, rhsI, bc] = apply_BC(P3, K, f, boundaryValue)
% [APPLY_BC]
%
% Apply essential Dirichlet boundary conditions to the P3 FEM system
%
%       K * psi = f
%
% by forming the reduced system
%
%       KII * psiI = fI - KIB * psiB.
%
%
% INPUT
%
%   P3
%       P3 mesh structure containing
%
%           P3.points         : Ndof x 2 coordinates [R,Z]
%           P3.boundaryNodes  : boundary-node indices
%
%   K
%       Ndof x Ndof global stiffness matrix
%
%   f
%       Ndof x 1 global source vector
%
%   boundaryValue
%       Prescribed physical poloidal flux on the boundary.
%
%       Supported forms:
%
%           scalar
%               Same value at every boundary node.
%
%           Nb x 1 vector
%               Values ordered according to P3.boundaryNodes.
%
%           Ndof x 1 vector
%               Full nodal vector; boundary entries are selected.
%
%           function handle
%               boundaryValue(R,Z), evaluated at boundary nodes.
%
%
% OUTPUT
%
%   KII
%       Reduced stiffness matrix on free nodes.
%
%   rhsI
%       Dirichlet-corrected right-hand side.
%
%   bc
%       Structure containing
%
%           bc.freeNodes
%           bc.boundaryNodes
%           bc.boundaryValues
%           bc.Ndof
    arguments
        P3 (1,1) struct
        K (:, :) double
        f (:, 1) double
        boundaryValue (1,1) double = 0.0
    end

    points = P3.points;
    Ndof = size(points, 1);

    %==============================================
    % Matrix and source-vector validation
    %==============================================
    if ~isnumeric(K) || ~isreal(K) || ~isequal(size(K), [Ndof,Ndof])
        error(['[APPLY_BC] K must be real ', 'Ndof x Ndof matrix.'])
    end

    if any(~isfinite(nonzeros(K)))
        error('[APPLY_BC] K contains NaN or Inf.');
    end

    if ~isnumeric(f) || ~isreal(f) || ~isvector(f) || numel(f) ~= Ndof
        error(['[APPLY_BC] f must contain', 'one real value per global DOF.']);
    end

    f = f(:);

    if any(~isfinite(f))
        error('[APPLY_BC] f contains NaN or Inf.');
    end

    %==============================================
    % Boundary-node validation
    %==============================================
    boundaryNodes = P3.boundaryNodes(:);

    if isempty(boundaryNodes) || any(~isfinite(boundaryNodes)) || ...
            any(boundaryNodes ~= round(boundaryNodes)) || ...
            any(boundaryNodes < 1) || ...
            any(boundaryNodes > Ndof)
        error('[APPLY_BC] Invalid boundary-node indices.');
    end

    if numel(unique(boundaryNodes)) ~= numel(boundaryNodes)
        error('[APPLY_BC] boundaryNodes contains duplicated indices.');
    end

    Nb = numel(boundaryNodes);

    %==============================================
    % Evaluate prescribed boundary flux
    %==============================================
    psiB = repmat(boundaryValue, Nb, 1);

    %==============================================
    % Identify free nodes
    %==============================================
    isBoundary = false(Ndof, 1);
    isBoundary(boundaryNodes) = true;

    freeNodes = find(~isBoundary);

    if isempty(freeNodes)
        error('[APPLY_BC] No free DOFs remain.');
    end

    %==============================================
    % Reduced Dirichlet system
    %==============================================
    KII = K(freeNodes, freeNodes);
    KIB = K(freeNodes, boundaryNodes);

    rhsI = f(freeNodes) - KIB*psiB;

    %==============================================
    % Information required to restore
    % the full solution
    %==============================================
    bc = struct();

    bc.Ndof           = Ndof;
    bc.freeNodes      = freeNodes;
    bc.boundaryNodes  = boundaryNodes;
    bc.boundaryValues = psiB;
end