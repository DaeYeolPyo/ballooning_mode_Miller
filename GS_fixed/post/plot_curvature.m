function fig = plot_curvature(meshData,curvature,opts)
%PLOT_CURVATURE Plot kappa dot grad(p) and its zero-level boundaries.

narginchk(2,3);
if nargin<3 || isempty(opts)
    opts = struct();
end
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidCurvaturePlotOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct('titlePrefix','\kappa\cdot\nabla p', ...
    'showMesh',false);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownCurvaturePlotOption', ...
        'Unknown plot option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
if isstring(opts.titlePrefix) && isscalar(opts.titlePrefix)
    opts.titlePrefix = char(opts.titlePrefix);
end
if ~ischar(opts.titlePrefix) || ~isrow(opts.titlePrefix)
    error('GS:post:InvalidCurvaturePlotTitle', ...
        'titlePrefix must be text.');
end
opts.showMesh = logical(opts.showMesh);

requiredMesh = {'nodes','elements','boundaryPoints'};
if ~isstruct(meshData) || ~all(isfield(meshData,requiredMesh)) || ...
        ~isstruct(curvature) || ~isfield(curvature,'node') || ...
        ~isfield(curvature.node,'kappaDotGradP')
    error('GS:post:InvalidCurvaturePlotInput', ...
        'Use a mesh and result returned by compute_magnetic_curvature.');
end

if opts.showMesh
    edgeColor = [0.75,0.75,0.75];
else
    edgeColor = 'none';
end

values = curvature.node.kappaDotGradP(:);
if numel(values)~=size(meshData.nodes,1) || any(~isfinite(values))
    error('GS:post:InvalidCurvaturePlotInput', ...
        'The nodal kappaDotGradP field is invalid.');
end
zeroSet = extract_flux_contours(meshData,values,0);

fig = figure('Name','Kappa dot grad(p)','Color','w');
ax = axes(fig);
patch(ax,'Faces',meshData.elements,'Vertices',meshData.nodes, ...
    'FaceVertexCData',values,'FaceColor','interp','EdgeColor',edgeColor);
hold(ax,'on')

boundary = meshData.boundaryPoints;
closedBoundary = boundary([1:end,1],:);
plot(ax,closedBoundary(:,1),closedBoundary(:,2), ...
    'k-','LineWidth',1.5,'DisplayName','LCFS')

zeroCurves = zeroSet.curves;
for k = 1:numel(zeroCurves)
    points = zeroCurves(k).points;
    % White underlay keeps the signed zero boundary visible over both ends
    % of the diverging color scale.
    plot(ax,points(:,1),points(:,2),'w-', ...
        'LineWidth',3.5,'HandleVisibility','off')
    if k==1
        visibility = 'on';
    else
        visibility = 'off';
    end
    plot(ax,points(:,1),points(:,2),'k--', ...
        'LineWidth',1.8,'HandleVisibility',visibility, ...
        'DisplayName','\kappa\cdot\nabla p = 0')
end
hold(ax,'off')

axis(ax,'equal','tight')
grid(ax,'on')
xlabel(ax,'R')
ylabel(ax,'Z')
title(ax,opts.titlePrefix)
colorbar(ax)
colormap(ax,blue_white_red(257))
limit = max(abs(values));
if limit>0
    clim(ax,[-limit,limit])
end
if ~isempty(zeroCurves)
    legend(ax,'Location','bestoutside')
end

fig.UserData = struct('zeroContours',zeroSet, ...
    'nZeroCurves',zeroSet.nCurves,'nClosedZeroCurves',zeroSet.nClosed, ...
    'nOpenZeroCurves',zeroSet.nOpen);
end


function map = blue_white_red(n)
half = floor(n/2);
blue = [0.10,0.25,0.85];
white = [1,1,1];
red = [0.85,0.15,0.10];
lower = [linspace(blue(1),white(1),half+1).', ...
         linspace(blue(2),white(2),half+1).', ...
         linspace(blue(3),white(3),half+1).'];
upperCount = n-size(lower,1);
upper = [linspace(white(1),red(1),upperCount+1).', ...
         linspace(white(2),red(2),upperCount+1).', ...
         linspace(white(3),red(3),upperCount+1).'];
map = [lower;upper(2:end,:)];
end
