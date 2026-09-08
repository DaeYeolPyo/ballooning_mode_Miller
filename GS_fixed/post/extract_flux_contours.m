function contourSets = extract_flux_contours(meshData,psi,levels,opts)
%EXTRACT_FLUX_CONTOURS Trace P1 level sets directly on triangular elements.
%
%   SETS = EXTRACT_FLUX_CONTOURS(MESH,PSI,LEVELS) intersects each triangle
%   with every requested flux level, then stitches the segments into ordered
%   open or closed polylines. This avoids interpolation across the concave
%   region outside a C-shaped plasma.
%
%   OPTS fields:
%       pointTolerance : endpoint stitching tolerance; default is 1e-9 of
%                        the mesh diameter

%   SETS(k) contains level, curves, nCurves, nClosed, nOpen, and nSegments.

narginchk(3,4);
[nodes,elements] = validate_inputs(meshData,psi);
psi = psi(:);
validateattributes(levels,{'numeric'}, ...
    {'real','finite','vector','nonempty'},mfilename,'levels');
levels = levels(:);

if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts,nodes);

emptyCurve = struct('points',{},'closed',{},'length',{}, ...
    'signedArea',{},'orientation',{});
contourSets = repmat(struct('level',[],'curves',emptyCurve, ...
    'nCurves',0,'nClosed',0,'nOpen',0,'nSegments',0),numel(levels),1);

elementValues = reshape(psi(elements(:)),size(elements));
elementScale = max(1,max(abs(psi)));
valueTolerance = 32*eps(elementScale);

for levelIndex = 1:numel(levels)
    level = levels(levelIndex);
    segmentStart = zeros(size(elements,1),2);
    segmentEnd = zeros(size(elements,1),2);
    nSegments = 0;

    for elementIndex = 1:size(elements,1)
        value = elementValues(elementIndex,:);
        if level < min(value)-valueTolerance || ...
                level > max(value)+valueTolerance
            continue
        end
        signed = value-level;
        signed(abs(signed)<=valueTolerance) = valueTolerance;
        localEdges = [1,2;2,3;3,1];
        intersections = zeros(3,2);
        nIntersections = 0;
        for edgeIndex = 1:3
            a = localEdges(edgeIndex,1);
            b = localEdges(edgeIndex,2);
            if signed(a)*signed(b) < 0
                denominator = value(b)-value(a);
                t = (level-value(a))/denominator;
                pointA = nodes(elements(elementIndex,a),:);
                pointB = nodes(elements(elementIndex,b),:);
                nIntersections = nIntersections+1;
                intersections(nIntersections,:) = pointA+t*(pointB-pointA);
            end
        end
        if nIntersections == 2
            nSegments = nSegments+1;
            segmentStart(nSegments,:) = intersections(1,:);
            segmentEnd(nSegments,:) = intersections(2,:);
        end
    end

    segmentStart = segmentStart(1:nSegments,:);
    segmentEnd = segmentEnd(1:nSegments,:);
    curves = stitch_segments(segmentStart,segmentEnd,opts.pointTolerance);
    closedFlags = [curves.closed];

    contourSets(levelIndex).level = level;
    contourSets(levelIndex).curves = curves;
    contourSets(levelIndex).nCurves = numel(curves);
    contourSets(levelIndex).nClosed = nnz(closedFlags);
    contourSets(levelIndex).nOpen = numel(curves)-nnz(closedFlags);
    contourSets(levelIndex).nSegments = nSegments;
end

end


function curves = stitch_segments(startPoint,endPoint,tolerance)
emptyCurve = struct('points',{},'closed',{},'length',{}, ...
    'signedArea',{},'orientation',{});
if isempty(startPoint)
    curves = emptyCurve;
    return
end

allPoints = [startPoint;endPoint];
keys = round(allPoints/tolerance);
[~,~,pointId] = unique(keys,'rows','stable');
nSegments = size(startPoint,1);
edgePointId = [pointId(1:nSegments),pointId(nSegments+1:end)];
valid = edgePointId(:,1) ~= edgePointId(:,2);
edgePointId = edgePointId(valid,:);
startPoint = startPoint(valid,:);
endPoint = endPoint(valid,:);
nSegments = size(edgePointId,1);
if nSegments == 0
    curves = emptyCurve;
    return
end

nPoints = max(edgePointId(:));
coordinate = zeros(nPoints,2);
count = zeros(nPoints,1);
for k = 1:nSegments
    ids = edgePointId(k,:);
    coordinate(ids(1),:) = coordinate(ids(1),:)+startPoint(k,:);
    coordinate(ids(2),:) = coordinate(ids(2),:)+endPoint(k,:);
    count(ids) = count(ids)+1;
end
coordinate = bsxfun(@rdivide,coordinate,count);

incident = cell(nPoints,1);
for k = 1:nSegments
    incident{edgePointId(k,1)}(end+1) = k;
    incident{edgePointId(k,2)}(end+1) = k;
end

used = false(nSegments,1);
curves = emptyCurve;
while any(~used)
    unusedEdges = find(~used);
    endpoints = unique(edgePointId(unusedEdges,:));
    degree = zeros(nPoints,1);
    for k = unusedEdges.'
        degree(edgePointId(k,:)) = degree(edgePointId(k,:))+1;
    end
    openCandidates = endpoints(degree(endpoints)==1);
    if isempty(openCandidates)
        current = min(endpoints);
    else
        current = min(openCandidates);
    end
    path = current;
    previousEdge = 0;
    while true
        candidates = incident{current};
        candidates = candidates(~used(candidates));
        if previousEdge ~= 0
            candidates(candidates==previousEdge) = [];
        end
        if isempty(candidates)
            break
        end
        edge = candidates(1);
        used(edge) = true;
        endpointsOfEdge = edgePointId(edge,:);
        next = endpointsOfEdge(endpointsOfEdge~=current);
        path(end+1,1) = next; %#ok<AGROW>
        previousEdge = edge;
        current = next;
        if current == path(1)
            break
        end
    end

    points = coordinate(path,:);
    closed = numel(path)>2 && path(end)==path(1);
    delta = diff(points,1,1);
    lengthValue = sum(sqrt(sum(delta.^2,2)));
    if closed
        areaValue = 0.5*sum(points(1:end-1,1).*points(2:end,2) ...
                            -points(2:end,1).*points(1:end-1,2));
        if areaValue >= 0
            orientation = 'counterclockwise';
        else
            orientation = 'clockwise';
        end
    else
        areaValue = NaN;
        orientation = 'open';
    end
    curve = struct('points',points,'closed',closed, ...
        'length',lengthValue,'signedArea',areaValue, ...
        'orientation',orientation);
    curves(end+1,1) = curve; %#ok<AGROW>
end
end


function [nodes,elements] = validate_inputs(meshData,psi)
required = {'nodes','elements'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData,required))
    error('GS:post:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
nodes = meshData.nodes;
elements = meshData.elements;
validateattributes(psi,{'numeric'}, ...
    {'real','finite','vector','numel',size(nodes,1)},mfilename,'psi');
end


function opts = validate_options(opts,nodes)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidContourOptions','opts must be a scalar struct.');
end
diameter = hypot(max(nodes(:,1))-min(nodes(:,1)), ...
                 max(nodes(:,2))-min(nodes(:,2)));
defaults = struct('pointTolerance',max(1e-12,1e-9*diameter));
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownContourOption', ...
        'Unknown contour option: %s',strjoin(unknown,', '));
end
if ~isfield(opts,'pointTolerance')
    opts.pointTolerance = defaults.pointTolerance;
end
validateattributes(opts.pointTolerance,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'opts.pointTolerance');
end
