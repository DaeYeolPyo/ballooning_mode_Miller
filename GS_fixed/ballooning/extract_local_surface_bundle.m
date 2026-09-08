function bundle = extract_local_surface_bundle( ...
    meshData,psi,axisData,psiNValues,opts)
%EXTRACT_LOCAL_SURFACE_BUNDLE Extract periodic FEM flux-surface geometry.
%
%   BUNDLE = EXTRACT_LOCAL_SURFACE_BUNDLE(MESH,PSI,AXIS,PSINVALUES)
%   traces the requested nested P1 flux contours directly on the triangular
%   plasma mesh. Each surface is started at its outboard-midplane crossing,
%   traversed upward, and resampled at uniform geometric arclength using
%
%       vartheta = 2*pi*s/L,       0 <= vartheta < 2*pi.
%
%   This vartheta is only a periodic geometric coordinate. It is not yet a
%   straight-field-line angle. Building the PEST angle is a later stage.
%
%   OPTS fields:
%       nPoloidal      points per surface (default 256)
%       boundaryValue  LCFS psi (boundary median by default)
%
%   Output surface fields include R, Z, points, vartheta, s, dl, perimeter,
%   area, raw contour points, midplane crossings, and a P1 interpolation
%   error showing that the resampled curve remains on its requested psi.

narginchk(4,5);
if nargin<5 || isempty(opts)
    opts = struct();
end
[nodes,elements,boundaryNodes] = validate_mesh(meshData);
validateattributes(psi,{'numeric'}, ...
    {'real','finite','vector','numel',size(nodes,1)},mfilename,'psi');
psi = psi(:);
validate_axis(axisData);
validateattributes(psiNValues,{'numeric'}, ...
    {'real','finite','vector','nonempty','>',0,'<',1}, ...
    mfilename,'psiNValues');
psiNValues = psiNValues(:);
if any(diff(psiNValues)<=0)
    error('GS:ballooning:NonIncreasingPsiN', ...
        'psiNValues must be strictly increasing.');
end
opts = validate_options(opts,psi,boundaryNodes);

if exist('extract_flux_contours','file')~=2
    error('GS:ballooning:MissingContourDependency', ...
        'Add the post folder containing extract_flux_contours to the path.');
end

psiAxis = axisData.psi;
psiBoundary = opts.boundaryValue;
span = psiBoundary-psiAxis;
if abs(span)<=100*eps(max([1,abs(psiAxis),abs(psiBoundary)]))
    error('GS:ballooning:CollapsedFluxRange', ...
        'Axis and boundary flux values must differ.');
end
levels = psiAxis+psiNValues*span;
contourSets = extract_flux_contours(meshData,psi,levels);

nSurface = numel(psiNValues);
emptySurface = struct('psiN',[],'psi',[],'vartheta',[], ...
    'R',[],'Z',[],'points',[],'s',[],'dl',[],'perimeter',[], ...
    'area',[],'signedArea',[],'orientation',[], ...
    'outboardMidplaneR',[],'inboardMidplaneR',[], ...
    'rawPoints',[],'rawContourIndex',[], ...
    'maxAbsFluxError',[],'rmsFluxError',[], ...
    'relativeChordVariation',[],'axisEnclosed',[]);
surfaces = repmat(emptySurface,nSurface,1);

for k = 1:nSurface
    [rawPoints,curveIndex] = choose_axis_contour( ...
        contourSets(k),axisData.R,axisData.Z);
    rawPoints = remove_duplicate_endpoint(rawPoints);
    rawPoints = orient_outboard_up(rawPoints,axisData.Z);
    [points,s,dl,perimeter] = resample_closed_polyline( ...
        rawPoints,axisData.Z,opts.nPoloidal);

    sampledPsi = interpolate_p1(nodes,elements,psi,points);
    fluxError = sampledPsi-levels(k);
    signedArea = polygon_signed_area(points);
    [outboardR,inboardR] = midplane_crossings(points,axisData.Z);

    surfaces(k).psiN = psiNValues(k);
    surfaces(k).psi = levels(k);
    surfaces(k).vartheta = 2*pi*s/perimeter;
    surfaces(k).R = points(:,1);
    surfaces(k).Z = points(:,2);
    surfaces(k).points = points;
    surfaces(k).s = s;
    surfaces(k).dl = dl;
    surfaces(k).perimeter = perimeter;
    surfaces(k).area = abs(signedArea);
    surfaces(k).signedArea = signedArea;
    surfaces(k).orientation = orientation_name(signedArea);
    surfaces(k).outboardMidplaneR = outboardR;
    surfaces(k).inboardMidplaneR = inboardR;
    surfaces(k).rawPoints = rawPoints;
    surfaces(k).rawContourIndex = curveIndex;
    surfaces(k).maxAbsFluxError = max(abs(fluxError));
    surfaces(k).rmsFluxError = sqrt(mean(fluxError.^2));
    surfaces(k).relativeChordVariation = std(dl)/mean(dl);
    surfaces(k).axisEnclosed = inpolygon( ...
        axisData.R,axisData.Z,points(:,1),points(:,2));
end

areas = [surfaces.area].';
if any(diff(areas)<=0)
    error('GS:ballooning:NonNestedSurfaceAreas', ...
        'Flux-surface area must increase strictly with psiN.');
end
if ~all([surfaces.axisEnclosed])
    error('GS:ballooning:AxisNotEnclosed', ...
        'Every selected surface must enclose the magnetic axis.');
end

bundle = struct();
bundle.psiAxis = psiAxis;
bundle.psiBoundary = psiBoundary;
bundle.psiN = psiNValues;
bundle.psi = levels;
bundle.vartheta = surfaces(1).vartheta;
bundle.nSurface = nSurface;
bundle.nPoloidal = opts.nPoloidal;
bundle.axis = axisData;
bundle.surfaces = surfaces;
bundle.options = opts;
bundle.diagnostics = struct( ...
    'maximumFluxError',max([surfaces.maxAbsFluxError]), ...
    'maximumRelativeChordVariation', ...
        max([surfaces.relativeChordVariation]), ...
    'areasStrictlyIncreasing',all(diff(areas)>0), ...
    'allEncloseAxis',all([surfaces.axisEnclosed]), ...
    'allStartAtAxisHeight',max(abs( ...
        arrayfun(@(surface) surface.Z(1),surfaces)-axisData.Z)));
end


function [points,index] = choose_axis_contour(set,axisR,axisZ)
curves = set.curves;
candidate = zeros(0,1);
area = zeros(0,1);
for j = 1:numel(curves)
    if ~curves(j).closed
        continue
    end
    points = curves(j).points;
    if inpolygon(axisR,axisZ,points(:,1),points(:,2))
        candidate(end+1,1) = j; %#ok<AGROW>
        area(end+1,1) = abs(curves(j).signedArea); %#ok<AGROW>
    end
end
if isempty(candidate)
    error('GS:ballooning:NoAxisEnclosingContour', ...
        'No closed contour enclosing the magnetic axis was found.');
end
[~,largest] = max(area);
index = candidate(largest);
points = curves(index).points;
end


function points = remove_duplicate_endpoint(points)
scale = max([1;range(points(:,1));range(points(:,2))]);
if size(points,1)>1 && norm(points(end,:)-points(1,:))<1e-11*scale
    points(end,:) = [];
end
end


function points = orient_outboard_up(points,axisZ)
% Force counterclockwise first, then insert the outboard midplane crossing.
if polygon_signed_area(points)<0
    points = flipud(points);
end
points = start_at_outboard_crossing(points,axisZ);
if size(points,1)>1 && points(2,2)<points(1,2)
    points = flipud(points);
    points = start_at_outboard_crossing(points,axisZ);
end
end


function ordered = start_at_outboard_crossing(points,axisZ)
n = size(points,1);
crossingR = zeros(n,1);
crossingEdge = zeros(n,1);
crossingT = zeros(n,1);
nCrossing = 0;
for i = 1:n
    j = mod(i,n)+1;
    z1 = points(i,2)-axisZ;
    z2 = points(j,2)-axisZ;
    if (z1<=0 && z2>0) || (z1>0 && z2<=0)
        t = -z1/(z2-z1);
        nCrossing = nCrossing+1;
        crossingR(nCrossing) = points(i,1) ...
            +t*(points(j,1)-points(i,1));
        crossingEdge(nCrossing) = i;
        crossingT(nCrossing) = t;
    end
end
if nCrossing==0
    error('GS:ballooning:NoMidplaneCrossing', ...
        'The flux contour does not cross the magnetic-axis height.');
end
[~,choice] = max(crossingR(1:nCrossing));
i = crossingEdge(choice);
j = mod(i,n)+1;
t = crossingT(choice);
start = points(i,:)+t*(points(j,:)-points(i,:));
if j==1
    order = 1:n;
else
    order = [j:n,1:i];
end
ordered = [start;points(order,:)];
ordered = remove_duplicate_endpoint(ordered);
end


function [points,s,dl,L] = resample_closed_polyline(raw,axisZ,nPoints)
closed = [raw;raw(1,:)];
segmentLength = sqrt(sum(diff(closed,1,1).^2,2));
cumulative = [0;cumsum(segmentLength)];
keep = [true;diff(cumulative)>10*eps(max(1,cumulative(end)))];
cumulative = cumulative(keep);
closed = closed(keep,:);
L = cumulative(end);
s = (0:nPoints-1).'*L/nPoints;
R = interp1(cumulative,closed(:,1),s,'linear');
Z = interp1(cumulative,closed(:,2),s,'linear');
points = [R,Z];

% Roundoff in interpolation can move the first Z by a few ulps.
points(1,2) = axisZ;
delta = circshift(points,-1,1)-points;
dl = sqrt(sum(delta.^2,2));
end


function values = interpolate_p1(nodes,elements,nodalValues,points)
tri = triangulation(elements,nodes);
triangleIndex = pointLocation(tri,points);
if any(isnan(triangleIndex))
    error('GS:ballooning:ResamplingOutsideMesh', ...
        ['A resampled contour point fell outside the FEM mesh. Increase ' ...
         'the contour resolution or inspect the selected surface.']);
end
barycentric = cartesianToBarycentric(tri,triangleIndex,points);
localValues = reshape( ...
    nodalValues(elements(triangleIndex,:)),size(barycentric));
values = sum(barycentric.*localValues,2);
end


function [outboardR,inboardR] = midplane_crossings(points,axisZ)
n = size(points,1);
values = zeros(n,1);
count = 0;
for i = 1:n
    j = mod(i,n)+1;
    z1 = points(i,2)-axisZ;
    z2 = points(j,2)-axisZ;
    if (z1<=0 && z2>0) || (z1>0 && z2<=0)
        t = -z1/(z2-z1);
        count = count+1;
        values(count) = points(i,1)+t*(points(j,1)-points(i,1));
    end
end
values = values(1:count);
outboardR = max(values);
inboardR = min(values);
end


function value = polygon_signed_area(points)
next = points([2:end,1],:);
value = 0.5*sum(points(:,1).*next(:,2)-next(:,1).*points(:,2));
end


function name = orientation_name(signedArea)
if signedArea>0
    name = 'counterclockwise';
else
    name = 'clockwise';
end
end


function [nodes,elements,boundaryNodes] = validate_mesh(meshData)
required = {'nodes','elements','boundaryNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData,required))
    error('GS:ballooning:InvalidMesh', ...
        'meshData must be returned by generate_mesh.');
end
nodes = meshData.nodes;
elements = meshData.elements;
boundaryNodes = meshData.boundaryNodes(:);
if size(nodes,2)~=2 || size(elements,2)~=3 || ...
        any(~isfinite(nodes(:))) || any(nodes(:,1)<=0)
    error('GS:ballooning:InvalidMesh','Mesh arrays are invalid.');
end
end


function validate_axis(axisData)
required = {'R','Z','psi'};
if ~isstruct(axisData) || ~isscalar(axisData) || ...
        ~all(isfield(axisData,required)) || ...
        any(~isfinite([axisData.R,axisData.Z,axisData.psi]))
    error('GS:ballooning:InvalidAxis', ...
        'axisData must be returned by find_axis.');
end
end


function opts = validate_options(opts,psi,boundaryNodes)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:ballooning:InvalidSurfaceOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('nPoloidal',256, ...
    'boundaryValue',median(psi(boundaryNodes)));
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:ballooning:UnknownSurfaceOption', ...
        'Unknown surface option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.nPoloidal,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',32}, ...
    mfilename,'opts.nPoloidal');
validateattributes(opts.boundaryValue,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.boundaryValue');
end
