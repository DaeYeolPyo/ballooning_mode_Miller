function mesh = generate_mesh(geom, opts)
%GENERATE_MESH Generate a constrained P1 triangular mesh inside an LCFS.
%
%   MESH = GENERATE_MESH(GEOM) meshes the simple polygon returned by
%   MAKE_BOUNDARY. A staggered triangular point lattice supplies interior
%   nodes, and the LCFS edges are imposed as Delaunay constraints.
%
%   MESH = GENERATE_MESH(GEOM, OPTS) supports:
%       targetH             nominal interior-node spacing
%       boundaryClearance   excluded distance from the LCFS, in units of h
%       maxGridPoints       safety limit before triangulation
%
%   The output contains nodes, P1 elements, LCFS boundary edges/nodes,
%   interior nodes, element areas, and basic mesh-quality diagnostics.

narginchk(1, 2);

required = {'points', 'area', 'characteristicLength', 'nBoundary'};
if ~isstruct(geom) || ~isscalar(geom) || ...
        ~all(isfield(geom, required))
    error('GS:geometry:InvalidGeometry', ...
        'geom must be a scalar struct returned by make_boundary.');
end

if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:geometry:InvalidMeshOptions', ...
        'opts must be a scalar struct.');
end

defaults = struct( ...
    'targetH', geom.characteristicLength/20, ...
    'boundaryClearance', 0.12, ...
    'maxGridPoints', 2e6);

unknown = setdiff(fieldnames(opts), fieldnames(defaults));
if ~isempty(unknown)
    error('GS:geometry:UnknownMeshOption', ...
        'Unknown mesh option: %s', strjoin(unknown, ', '));
end

names = fieldnames(defaults);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(opts, name)
        opts.(name) = defaults.(name);
    end
end

validateattributes(opts.targetH, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'opts.targetH');
validateattributes(opts.boundaryClearance, {'numeric'}, ...
    {'real', 'finite', 'scalar', '>=', 0, '<=', 0.5}, ...
    mfilename, 'opts.boundaryClearance');
validateattributes(opts.maxGridPoints, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
    mfilename, 'opts.maxGridPoints');

boundary = geom.points;
nBoundary = size(boundary, 1);
if nBoundary ~= geom.nBoundary || size(boundary, 2) ~= 2
    error('GS:geometry:InconsistentGeometry', ...
        'geom.points and geom.nBoundary are inconsistent.');
end

h = opts.targetH;
dy = sqrt(3)*h/2;
xMin = min(boundary(:,1));
xMax = max(boundary(:,1));
yMin = min(boundary(:,2));
yMax = max(boundary(:,2));

estimatedGridPoints = ceil((xMax-xMin)/h + 3) ...
                    * ceil((yMax-yMin)/dy + 3);
if estimatedGridPoints > opts.maxGridPoints
    error('GS:geometry:MeshTooLarge', ...
        ['The requested targetH would create approximately %.0f candidate ' ...
         'grid points, above maxGridPoints = %.0f.'], ...
        estimatedGridPoints, opts.maxGridPoints);
end

yRows = (yMin-dy:dy:yMax+dy).';
rowPoints = cell(numel(yRows), 1);
for row = 1:numel(yRows)
    offset = 0.5*h*mod(row-1, 2);
    xRow = (xMin-h+offset:h:xMax+h).';
    yRow = repmat(yRows(row), size(xRow));
    [inside, onBoundary] = inpolygon( ...
        xRow, yRow, boundary(:,1), boundary(:,2));
    keep = inside & ~onBoundary;
    rowPoints{row} = [xRow(keep), yRow(keep)];
end
interior = vertcat(rowPoints{:});

if ~isempty(interior) && opts.boundaryClearance > 0
    distance = distance_to_boundary(interior, boundary);
    interior(distance < opts.boundaryClearance*h, :) = [];
end

if size(interior, 1) < 3
    error('GS:geometry:InsufficientInteriorPoints', ...
        ['targetH is too large for this LCFS. Reduce targetH until at least ' ...
         'three interior nodes are generated.']);
end

constraints = [(1:nBoundary).', [2:nBoundary, 1].'];
allPoints = [boundary; interior];
dt = delaunayTriangulation(allPoints, constraints);

nodes = dt.Points;
elements = dt.ConnectivityList;
centroids = (nodes(elements(:,1),:) + nodes(elements(:,2),:) ...
           + nodes(elements(:,3),:))/3;
[inside, onBoundary] = inpolygon(centroids(:,1), centroids(:,2), ...
    boundary(:,1), boundary(:,2));
elements = elements(inside | onBoundary, :);

twiceArea = signed_twice_area(nodes, elements);
areaTol = 100*eps(max(1, geom.area));
elements(abs(twiceArea) <= areaTol, :) = [];
twiceArea(abs(twiceArea) <= areaTol) = [];

clockwise = twiceArea < 0;
elements(clockwise, [2, 3]) = elements(clockwise, [3, 2]);

used = unique(elements(:));
oldToNew = zeros(size(nodes,1), 1);
oldToNew(used) = (1:numel(used)).';
nodes = nodes(used, :);
elements = oldToNew(elements);

elementArea = 0.5*signed_twice_area(nodes, elements);
if any(elementArea <= 0)
    error('GS:geometry:InvalidElementOrientation', ...
        'The generated mesh contains a degenerate or clockwise element.');
end

[boundaryEdges, boundaryOrder] = extract_boundary(elements);
boundaryNodes = boundaryOrder(:);
interiorNodes = setdiff((1:size(nodes,1)).', boundaryNodes, 'stable');

if numel(boundaryNodes) ~= nBoundary
    error('GS:geometry:BoundaryNotPreserved', ...
        ['The constrained triangulation did not preserve every sampled LCFS ' ...
         'node. Reduce targetH or inspect the boundary sampling.']);
end

meshArea = sum(elementArea);
relativeAreaError = abs(meshArea - geom.area)/geom.area;
if relativeAreaError > 1e-8
    error('GS:geometry:AreaMismatch', ...
        ['Triangulation area differs from polygon area by %.3e. This usually ' ...
         'indicates a disconnected or incorrectly filtered mesh.'], ...
        relativeAreaError);
end

quality = triangle_quality(nodes, elements, elementArea);

mesh = struct();
mesh.nodes = nodes;
mesh.elements = elements;
mesh.boundaryEdges = boundaryEdges;
mesh.boundaryNodes = boundaryNodes;
mesh.boundaryOrder = boundaryOrder;
mesh.boundaryPoints = nodes(boundaryOrder, :);
mesh.interiorNodes = interiorNodes;
mesh.elementArea = elementArea;
mesh.nNodes = size(nodes, 1);
mesh.nElements = size(elements, 1);
mesh.nBoundaryNodes = numel(boundaryNodes);
mesh.targetH = h;
mesh.options = opts;
mesh.area = meshArea;
mesh.quality = struct( ...
    'perElement', quality, ...
    'minimum', min(quality), ...
    'mean', mean(quality), ...
    'relativeAreaError', relativeAreaError);

end


function distance = distance_to_boundary(points, boundary)
nPoints = size(points, 1);
nBoundary = size(boundary, 1);
distanceSquared = inf(nPoints, 1);

for k = 1:nBoundary
    kNext = mod(k, nBoundary) + 1;
    a = boundary(k,:);
    edge = boundary(kNext,:) - a;
    edgeSquared = dot(edge, edge);
    if edgeSquared == 0
        continue
    end

    relative = points - a;
    projection = (relative*edge.')/edgeSquared;
    projection = max(0, min(1, projection));
    closest = a + projection.*edge;
    candidate = sum((points - closest).^2, 2);
    distanceSquared = min(distanceSquared, candidate);
end

distance = sqrt(distanceSquared);
end


function twiceArea = signed_twice_area(nodes, elements)
p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
twiceArea = (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) ...
          - (p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1));
end


function [boundaryEdges, order] = extract_boundary(elements)
allEdges = sort([elements(:,[1,2]); ...
                 elements(:,[2,3]); ...
                 elements(:,[3,1])], 2);
[uniqueEdges, ~, edgeMap] = unique(allEdges, 'rows');
edgeCount = accumarray(edgeMap, 1);
boundaryEdges = uniqueEdges(edgeCount == 1, :);

degree = accumarray(boundaryEdges(:), 1);
boundaryNodes = unique(boundaryEdges(:));
if any(degree(boundaryNodes) ~= 2)
    error('GS:geometry:InvalidBoundaryGraph', ...
        'The triangulation boundary is not a single closed polygonal graph.');
end

nEdges = size(boundaryEdges, 1);
adjacency = cell(max(boundaryNodes), 1);
for k = 1:nEdges
    i = boundaryEdges(k,1);
    j = boundaryEdges(k,2);
    adjacency{i}(end+1) = j;
    adjacency{j}(end+1) = i;
end

order = zeros(nEdges, 1);
order(1) = min(boundaryNodes);
previous = 0;
current = order(1);
for k = 2:nEdges
    candidates = adjacency{current};
    candidates(candidates == previous) = [];
    if isempty(candidates)
        error('GS:geometry:OpenBoundaryGraph', ...
            'Could not traverse the triangulation boundary as a closed loop.');
    end
    next = candidates(1);
    order(k) = next;
    previous = current;
    current = next;
end

if ~ismember(order(1), adjacency{order(end)}) || ...
        numel(unique(order)) ~= nEdges
    error('GS:geometry:MultipleBoundaryLoops', ...
        'The triangulation boundary contains multiple loops or a branch.');
end
end


function quality = triangle_quality(nodes, elements, elementArea)
p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
sumEdgeSquared = sum((p2-p1).^2, 2) ...
               + sum((p3-p2).^2, 2) ...
               + sum((p1-p3).^2, 2);
quality = 4*sqrt(3)*elementArea./sumEdgeSquared;
end
