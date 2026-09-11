function meshP3 = generate_P3_nodes(meshP1)
% [GENERATE_P3_NODES]
% Upgrade a linear triangular mesh to P3 FEM connectivity.
%
% INPUT
%   meshP1 : triangulation object
%
% OUTPUT
%   meshP3 : structure containing
%
%       .points
%           Np3 x 2 array of [R,Z] coordinates of all P3 nodes.
%
%       .elements
%           Nc x 10 P3 element connectivity.
%
%           Local node ordering:
%
%                    3
%                    *
%                   / \
%                  8   7
%                 /     \
%                9  10   6
%               /    *    \
%              *---*---*---*
%              1   4   5   2
%
%           Nodes 4,5 : edge 1 -> 2
%           Nodes 6,7 : edge 2 -> 3
%           Nodes 8,9 : edge 3 -> 1
%           Node 10   : interior
%
%       .vertices
%           Original P1 node indices.
%
%       .edges
%           Ne x 2 canonical unique-edge vertex connectivity.
%
%       .edgeNodes
%           Ne x 2 global P3 node numbers.
%           Orientation follows meshP3.edges(:,1) -> meshP3.edges(:,2).
%
%       .interiorNodes
%           Nc x 1 global node numbers for element interior nodes.
%
%       .boundaryEdges
%           Boundary vertex edges returned by freeBoundary(DT).
%
%       .boundaryNodes
%           All P3 nodes lying on the boundary.
%
%       .elementEdges
%           Nt x 3 global edge numbers corresponding to local edges
%           [1-2], [2-3], [3-1].
%
% NOTES
%   Geometry remains straight-sided / affine.
%   The additional P3 nodes are solution interpolation nodes.
%
%   For a genuinely curved P3 isoparametric boundary, boundary edge-node
%   coordinates should later be projected onto the actual boundary.

    arguments
        meshP1 triangulation
    end

    % Original P1 mesh (only triangle vertices)
    vert = meshP1.Points;
    conn = meshP1.ConnectivityList;

    Nv = size(vert, 1);
    Nc = size(conn, 1);

    % Ensure counter-clockwise element orientation

    for e = 1:Nc
        v = conn(e, :);

        x1 = vert(v(1), :);
        x2 = vert(v(2), :);
        x3 = vert(v(3), :);

        detJ = (x2(1) - x1(1))*(x3(2) - x1(2)) ...
            - (x3(1) - x1(1))*(x2(2) - x1(2));

        if detJ < 0
            conn(e, [2 3]) = conn(e, [3 2]);
        elseif detJ == 0
            error('Degenerate triangle detected at element %d', e);
        end
    end

    % Construct local edges
    edge12 = conn(:, [1 2]);
    edge23 = conn(:, [2 3]);
    edge31 = conn(:, [3 1]);

    edge = [edge12; edge23; edge31];

    edge_sorted = sort(edge, 2);

    [edge_unique, ~, edgeID] = unique(edge_sorted, 'rows', 'stable');

    Ne = size(edge_unique, 1);

    edgeID12 = edgeID(1:Nc);
    edgeID23 = edgeID(Nc+1:2*Nc);
    edgeID31 = edgeID(2*Nc+1:3*Nc);

    elementEdges = [edgeID12, edgeID23, edgeID31];

    %===============================================
    % Allocate P3 global node numbering
    %
    % Original vertices:
    %    1 ... Nv
    %
    % First edge node:
    %    Nv+1 ... Nv+Ne
    %
    % Second edge node:
    %    Nv+Ne+1 ... Nv+2Ne
    %
    % Interior nodes:
    %    Nv+2Ne+1 ...
    %===============================================

    edgeNode1 = Nv + (1:Ne)';
    edgeNode2 = Nv + Ne + (1:Ne)';
    edgeNodes = [edgeNode1, edgeNode2];

    interiorNodes = Nv + 2*Ne + (1:Nc)';

    Np3 = Nv + 2*Ne + Nc;

    P3 = zeros(Np3, 2);
    P3(1:Nv, :) = vert; % Original P1 nodes

    %===============================================
    % Generate P3 nodes on every unique edge
    % 
    % Edge orientation:
    %
    %    i ----------------- j
    %
    %    i ----a-------b---- j
    %
    %    a = (2/3) xi + (1/3) xj
    %    b = (1/3) xi + (2/3) xj
    %===============================================

    for k = 1:Ne
        i = edge_unique(k, 1);
        j = edge_unique(k, 2);

        xi = vert(i, :);
        xj = vert(j, :);

        P3(edgeNode1(k), :) = (2/3)*xi + (1/3)*xj;
        P3(edgeNode2(k), :) = (1/3)*xi + (2/3)*xj;
    end

    %===============================================
    % Generate interior node
    %===============================================

    for e = 1:Nc
        vert_ind = conn(e, :);

        P3(interiorNodes(e), :) = ...
            (vert(vert_ind(1), :) + vert(vert_ind(2), :) + vert(vert_ind(3), :))/3;
    end

    %===============================================
    % Construct total connectivity
    %===============================================

    conn_tot = zeros(Nc, 10);

    for e = 1:Nc
        i1 = conn(e, 1);
        i2 = conn(e, 2);
        i3 = conn(e, 3);

        % Vertices
        conn_tot(e, 1:3) = [i1, i2, i3];

        %===========================
        % Edge 1 -> 2
        %===========================

        k = edgeID12(e);

        if i1 < i2
            conn_tot(e, 4:5) = [edgeNode1(k), edgeNode2(k)];
        else
            conn_tot(e, 4:5) = [edgeNode2(k), edgeNode1(k)];
        end

        %===========================
        % Edge 2 -> 3
        %===========================

        k = edgeID23(e);

        if i2 < i3
            conn_tot(e, 6:7) = [edgeNode1(k), edgeNode2(k)];
        else
            conn_tot(e, 6:7) = [edgeNode2(k), edgeNode1(k)];
        end

        %===========================
        % Edge 3 -> 1
        %===========================

        k = edgeID31(e);

        if i3 < i1
            conn_tot(e, 8:9) = [edgeNode1(k), edgeNode2(k)];
        else
            conn_tot(e, 8:9) = [edgeNode2(k), edgeNode1(k)];
        end

        %===========================
        % Interior node
        %===========================

        conn_tot(e, 10) = interiorNodes(e);
    end

    %===============================================
    % Find boundary edges
    %===============================================

    FB = freeBoundary(meshP1);

    % Convert free-boundary edge orientation into canonical orientation
    FBsorted = sort(FB, 2);

    % Find corresponding global edge IDs.
    [tf, boundaryEdgeID] = ismember(FBsorted, edge_unique, 'rows');

    if any(~tf)
        error('Failed to map some freeBoundary edges to global edge list.');
    end

    %===============================================
    % Construct complete P3 boundary-node list
    %===============================================

    boundaryVertexNodes = unique(FB(:));

    boundaryEdgeNode1 = edgeNode1(boundaryEdgeID);
    boundaryEdgeNode2 = edgeNode2(boundaryEdgeID);

    boundaryNodes = unique([ ...
        boundaryVertexNodes; boundaryEdgeNode1; boundaryEdgeNode2 ...
        ]);

    %===============================================
    % Package outputs
    %===============================================

    meshP3.points        = P3;
    meshP3.elements      = conn_tot;

    meshP3.vertices      = (1:Nv)';
    meshP3.edges         = edge_unique;
    meshP3.edgeNodes     = edgeNodes;
    meshP3.interiorNodes = interiorNodes;

    meshP3.elementEdges  = elementEdges;

    meshP3.boundaryEdges = FB;
    meshP3.boundaryNodes = boundaryNodes;

    meshP3.P1points      = vert;
    meshP3.P1elements    = conn;
end