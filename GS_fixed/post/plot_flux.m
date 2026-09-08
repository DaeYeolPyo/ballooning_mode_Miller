function output = plot_flux(meshData,psi,opts)
%PLOT_FLUX Plot a P1 flux field and contours without filling outside the LCFS.
%
%   OUTPUT = PLOT_FLUX(MESH,PSI) draws the triangular solution, normalized
%   contours, LCFS, and refined magnetic axis. OPTS supports:
%       normalizedLevels : default 0.1:0.1:0.9
%       showMesh         : draw triangle edges (default false)
%       showAxis         : mark magnetic axis (default true)
%       axisData         : precomputed find_axis result (default [])
%       boundaryValue    : LCFS flux (boundary median by default)
%       parent           : target axes handle (default new figure)
%       title            : plot title
%
%   OUTPUT contains graphics handles, axis data, and extracted contours.

narginchk(2,3);
if nargin < 3 || isempty(opts)
    opts = struct();
end
required = {'nodes','elements','boundaryNodes','boundaryOrder'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData,required))
    error('GS:post:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
validateattributes(psi,{'numeric'}, ...
    {'real','finite','vector','numel',size(meshData.nodes,1)},mfilename,'psi');
psi = psi(:);
opts = validate_options(opts,median(psi(meshData.boundaryNodes)));

if isempty(opts.axisData)
    axisData = find_axis(meshData,psi, ...
        struct('boundaryValue',opts.boundaryValue));
else
    axisData = opts.axisData;
end
levels = axisData.psi+opts.normalizedLevels ...
    *(opts.boundaryValue-axisData.psi);
contours = extract_flux_contours(meshData,psi,levels);

if isempty(opts.parent)
    figureHandle = figure('Name','Fixed-boundary Grad-Shafranov flux', ...
        'Color','w');
    axesHandle = axes('Parent',figureHandle);
else
    axesHandle = opts.parent;
    figureHandle = ancestor(axesHandle,'figure');
end

if opts.showMesh
    edgeColor = [0.70,0.74,0.78];
else
    edgeColor = 'none';
end
surfaceHandle = trisurf(meshData.elements,meshData.nodes(:,1), ...
    meshData.nodes(:,2),psi,'Parent',axesHandle, ...
    'EdgeColor',edgeColor,'FaceColor','interp');
view(axesHandle,2)
hold(axesHandle,'on')

colorMap = parula(numel(levels));
contourHandles = gobjects(0);
for k = 1:numel(contours)
    for j = 1:numel(contours(k).curves)
        points = contours(k).curves(j).points;
        contourHandles(end+1,1) = plot(axesHandle, ...
            points(:,1),points(:,2),'-','Color',colorMap(k,:), ...
            'LineWidth',1.35); %#ok<AGROW>
    end
end

boundary = meshData.nodes(meshData.boundaryOrder,:);
boundary = [boundary;boundary(1,:)];
boundaryHandle = plot(axesHandle,boundary(:,1),boundary(:,2), ...
    'k-','LineWidth',2.0);
if opts.showAxis
    axisHandle = plot(axesHandle,axisData.R,axisData.Z,'rp', ...
        'MarkerFaceColor','r','MarkerSize',10);
else
    axisHandle = gobjects(0);
end
hold(axesHandle,'off')
axis(axesHandle,'equal')
axis(axesHandle,'tight')
grid(axesHandle,'on')
xlabel(axesHandle,'R')
ylabel(axesHandle,'Z')
title(axesHandle,opts.title)
colorbar(axesHandle)

output = struct();
output.figure = figureHandle;
output.axes = axesHandle;
output.surface = surfaceHandle;
output.contourLines = contourHandles;
output.boundary = boundaryHandle;
output.axisMarker = axisHandle;
output.axisData = axisData;
output.contours = contours;
output.levels = levels;

end


function opts = validate_options(opts,boundaryDefault)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidPlotOptions','opts must be a scalar struct.');
end
defaults = struct('normalizedLevels',(0.1:0.1:0.9).', ...
    'showMesh',false,'showAxis',true,'axisData',[], ...
    'boundaryValue',boundaryDefault,'parent',[], ...
    'title','Fixed-boundary flux surfaces');
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownPlotOption', ...
        'Unknown plot option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.normalizedLevels,{'numeric'}, ...
    {'real','finite','vector','nonempty','>',0,'<',1}, ...
    mfilename,'opts.normalizedLevels');
opts.normalizedLevels = opts.normalizedLevels(:);
validateattributes(opts.boundaryValue,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.boundaryValue');
opts.showMesh = logical_scalar(opts.showMesh,'showMesh');
opts.showAxis = logical_scalar(opts.showAxis,'showAxis');
if isstring(opts.title) && isscalar(opts.title)
    opts.title = char(opts.title);
end
if ~ischar(opts.title)
    error('GS:post:InvalidPlotTitle','opts.title must be text.');
end
end


function value = logical_scalar(value,name)
validateattributes(value,{'logical','numeric'}, ...
    {'real','finite','scalar'},mfilename,['opts.' name]);
value = logical(value);
end
