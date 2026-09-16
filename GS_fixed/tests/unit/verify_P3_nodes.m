clc
clear
close all

%%
test_delaunay;

P1 = mesh.Points;
T1 = mesh.ConnectivityList;

Nv = size(P1,1);
Nt = size(T1,1);

tol = 1.e-12;

%% Run P3 node generation
meshP3 = generate_P3_nodes(mesh);

P3 = meshP3.points;
T3 = meshP3.elements;

%% 1. Count unique P1 edges

Eall = [T1(:, [1 2]); T1(:, [2 3]); T1(:, [3 1])];

Eall = sort(Eall,2);
E = unique(Eall,'rows');

Ne = size(E,1);

expectedNp3 = Nv + 2*Ne + Nt;
actualNp3   = size(P3,1);

fprintf('\n========================================\n');
fprintf(' P1 -> P3 mesh validation\n');
fprintf('========================================\n');

fprintf('Number of vertices       : %d\n', Nv);
fprintf('Number of edges          : %d\n', Ne);
fprintf('Number of triangles      : %d\n', Nt);
fprintf('Expected P3 nodes        : %d\n', expectedNp3);
fprintf('Actual P3 nodes          : %d\n', actualNp3);

assert(actualNp3 == expectedNp3, ...
    'FAIL: Wrong total number of P3 nodes.');

assert(size(T3,1) == Nt, ...
    'FAIL: Number of P3 elements differs from P1.');

assert(size(T3,2) == 10, ...
    'FAIL: A P3 triangle must contain 10 local nodes.');

fprintf('[PASS] Global node/element counts\n');

%% 2. Check node indices
assert(all(T3(:) >= 1), ...
    'FAIL: Invalid node index <= 0.');

assert(all(T3(:) <= actualNp3), ...
    'FAIL: Connectivity points outside P3 node array.');

assert(all(T3(:) == round(T3(:))), ...
    'FAIL: Connectivity contains non-integer indices.');

% No duplicated local node inside one element
for e = 1:Nt
    assert(numel(unique(T3(e,:))) == 10, ...
        'FAIL: Element %d has duplicated local nodes.', e);
end

fprintf('[PASS] Connectivity indices\n');

%% 3. Original P1 vertices must remain unchanged
vertexError = max(abs(P3(1:Nv,:) - P1), [], 'all');

fprintf('Maximum original-vertex error = %.3e\n', vertexError);

assert(vertexError < tol, ...
    'FAIL: Original P1 vertices were modified.');

fprintf('[PASS] Original P1 vertices preserved\n');

%% 4. First 3 local nodes must correspond to same P1 triangle

for e = 1:Nt

    originalVertices = sort(T1(e,:));
    p3Vertices       = sort(T3(e,1:3));

    assert(isequal(originalVertices, p3Vertices), ...
        'FAIL: Element %d has wrong vertex connectivity.', e);

end

fprintf('[PASS] P1/P3 element vertex correspondence\n');

%% 5. Check CCW orientation
minArea2 = inf;

for e = 1:Nt

    ids = T3(e,1:3);

    x1 = P3(ids(1),:);
    x2 = P3(ids(2),:);
    x3 = P3(ids(3),:);

    area2 = ...
        (x2(1)-x1(1))*(x3(2)-x1(2)) ...
      - (x2(2)-x1(2))*(x3(1)-x1(1));

    minArea2 = min(minArea2,area2);

    assert(area2 > tol, ...
        'FAIL: Element %d is not CCW or is degenerate.',e);

end

fprintf('Minimum 2*triangle area = %.3e\n',minArea2);
fprintf('[PASS] All elements CCW\n');

%% 6. Check geometric locations of all 10 local P3 nodes
maxGeometryError = 0;

for e = 1:Nt

    ids = T3(e,:);

    v1 = P3(ids(1),:);
    v2 = P3(ids(2),:);
    v3 = P3(ids(3),:);

    expected = [
        v1
        v2
        v3
        (2*v1 + v2)/3
        (v1 + 2*v2)/3
        (2*v2 + v3)/3
        (v2 + 2*v3)/3
        (2*v3 + v1)/3
        (v3 + 2*v1)/3
        (v1 + v2 + v3)/3
    ];

    actual = P3(ids,:);

    err = max(abs(actual-expected), [], 'all');

    maxGeometryError = max(maxGeometryError,err);

    assert(err < tol, ...
        'FAIL: Wrong P3 node position in element %d.',e);

end

fprintf('Maximum local geometry error = %.3e\n', ...
    maxGeometryError);

fprintf('[PASS] P3 local-node geometry\n');

%% 7. Barycentric-coordinate test
TRcheck = triangulation(T3(:,1:3),P3);

lambdaExpected = [
    1   0   0
    0   1   0
    0   0   1
    2/3 1/3 0
    1/3 2/3 0
    0   2/3 1/3
    0   1/3 2/3
    1/3 0   2/3
    2/3 0   1/3
    1/3 1/3 1/3
];

maxBaryError = 0;

for e = 1:Nt

    ids = T3(e,:);

    Xlocal = P3(ids,:);

    lambda = cartesianToBarycentric( ...
        TRcheck, ...
        repmat(e,10,1), ...
        Xlocal);

    err = max(abs(lambda-lambdaExpected),[],'all');

    maxBaryError = max(maxBaryError,err);

    assert(err < 10*tol, ...
        'FAIL: Barycentric coordinates wrong in element %d.',e);

end

fprintf('Maximum barycentric error = %.3e\n',maxBaryError);
fprintf('[PASS] Barycentric coordinates\n');

%% 8. Check sharing of P3 edge nodes
edgeMap = containers.Map( ...
    'KeyType','char', ...
    'ValueType','any');

for e = 1:Nt

    t = T3(e,:);

    % [edge vertex 1, edge vertex 2, P3 node a, P3 node b]
    localEdges = [
        t(1), t(2), t(4), t(5)
        t(2), t(3), t(6), t(7)
        t(3), t(1), t(8), t(9)
    ];

    for k = 1:3

        vertexPair = sort(localEdges(k,1:2));
        edgeNodes  = sort(localEdges(k,3:4));

        key = sprintf('%d_%d', ...
            vertexPair(1),vertexPair(2));

        if isKey(edgeMap,key)

            oldEdgeNodes = edgeMap(key);

            assert(isequal(oldEdgeNodes,edgeNodes), ...
                ['FAIL: Adjacent elements do not share ', ...
                 'the same P3 edge nodes on edge (%d,%d).'], ...
                 vertexPair(1),vertexPair(2));

        else

            edgeMap(key) = edgeNodes;

        end
    end
end

assert(edgeMap.Count == Ne, ...
    'FAIL: Number of reconstructed P3 edges is incorrect.');

fprintf('[PASS] Shared-edge global node consistency\n');

%% 9. Interior nodes must be unique to each triangle
interiorNodes = T3(:,10);

assert(numel(unique(interiorNodes)) == Nt, ...
    'FAIL: Interior node is shared by multiple triangles.');

fprintf('[PASS] Interior-node uniqueness\n');

%% 10. Every generated global node should actually be referenced
usedNodes = unique(T3(:));

assert(numel(usedNodes) == actualNp3, ...
    'FAIL: Some generated P3 nodes are never used.');

fprintf('[PASS] Every global P3 node is referenced\n');

%% 11. Connectivity visualization
% Settings
seedElement = ceil(Nt*rand());
nRing = 5;

% Build P1 triangulation from first 3 nodes of T3
TR = triangulation(T3(:, 1:3), P3);

% Find n-ring neighboring elements
selected = seedElement;
frontier = seedElement;

for ring = 1:nRing

    N = neighbors(TR,frontier);

    N = N(:);
    N = N(~isnan(N));

    N = unique(N);

    newElements = setdiff(N,selected);

    selected = unique([selected(:); newElements(:)]);
    frontier = newElements;

    if isempty(frontier)
        break
    end

end

% Draw selected triangle boundaries
figure;
hold on;
box on;

xlabel('R');
ylabel('Z');

title(sprintf( ...
    'P3 connectivity patch: seed E%d, %d-ring', ...
    seedElement,nRing));

for e = selected'

    ids = T3(e,1:3);

    X = P3(ids,:);

    Xclosed = [
        X
        X(1,:)
    ];

    if e == seedElement

        patch( ...
            X(:,1),X(:,2), ...
            [0.9 0.9 0.9], ...
            'FaceAlpha',0.5, ...
            'EdgeColor','k', ...
            'LineWidth',2);

    else

        plot( ...
            Xclosed(:,1), ...
            Xclosed(:,2), ...
            'k-', ...
            'LineWidth',1);

    end

end

% Collect nodes appearing in selected patch
patchNodes = unique(T3(selected,:));


vertexNodes = intersect( ...
    patchNodes, ...
    unique(T3(selected,1:3)));

edgeNodes = intersect( ...
    patchNodes, ...
    unique(T3(selected,4:9)));

interiorNodes = intersect( ...
    patchNodes, ...
    unique(T3(selected,10)));

scatter( ...
    P3(vertexNodes,1), ...
    P3(vertexNodes,2), ...
    80,'o','filled');

scatter( ...
    P3(edgeNodes,1), ...
    P3(edgeNodes,2), ...
    60,'s','filled');

scatter( ...
    P3(interiorNodes,1), ...
    P3(interiorNodes,2), ...
    70,'^','filled');

% Global node labels
xrange = max(P3(patchNodes,1)) - min(P3(patchNodes,1));
yrange = max(P3(patchNodes,2)) - min(P3(patchNodes,2));

scale = max(xrange,yrange);

if scale == 0
    scale = 1;
end

for g = patchNodes'

    text( ...
        P3(g,1) + 0.01*scale, ...
        P3(g,2) + 0.01*scale, ...
        sprintf('G%d',g), ...
        'FontSize',5, ...
        'FontWeight','bold');

end

% Element numbers
for e = selected'

    ids = T3(e,1:3);

    centroid = mean(P3(ids,:),1);

    text( ...
        centroid(1), ...
        centroid(2), ...
        sprintf('E%d',e), ...
        'HorizontalAlignment','center', ...
        'FontSize',7, ...
        'FontWeight','bold');

end