function axisData = find_axis(meshData, psi, opts)
%FIND_AXIS Locate and locally refine the magnetic-axis flux extremum.
%
%   AXISDATA = FIND_AXIS(MESHDATA, PSI) finds the interior nodal extremum
%   farthest in flux from the LCFS value, then fits a local quadratic to
%   refine (R_axis,Z_axis). The nodal result is retained if the fit is not a
%   definite extremum inside the plasma.
%
%   OPTS fields:
%       mode          'auto', 'min', or 'max' (default 'auto')
%       boundaryValue scalar LCFS flux; default is the boundary-node median
%       fitRings      neighboring graph rings used in the fit (default 2)

%   AXISDATA reports both nodal and refined estimates, the fitted Hessian,
%   gradient, classification, and whether refinement was accepted.

narginchk(2, 3);
[nodes, elements, boundaryNodes, interiorNodes] = validate_inputs( ...
    meshData, psi);
psi = psi(:);

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts, psi(boundaryNodes));

interiorValues = psi(interiorNodes);
[minimum, minimumIndex] = min(interiorValues);
[maximum, maximumIndex] = max(interiorValues);
mode = choose_mode(opts.mode, minimum, maximum, opts.boundaryValue);
if strcmp(mode, 'min')
    nodeValue = minimum;
    node = interiorNodes(minimumIndex);
else
    nodeValue = maximum;
    node = interiorNodes(maximumIndex);
end
nodePosition = nodes(node,:);

neighbors = vertex_neighbors(elements, size(nodes,1));
patch = graph_patch(node, neighbors, opts.fitRings);
fitAccepted = false;
fitPosition = nodePosition;
fitValue = nodeValue;
gradient = [NaN; NaN];
hessian = NaN(2);
fitResidual = NaN;

if numel(patch) >= 6
    displacement = nodes(patch,:)-nodePosition;
    A = [ones(numel(patch),1), displacement(:,1), displacement(:,2), ...
         0.5*displacement(:,1).^2, ...
         displacement(:,1).*displacement(:,2), ...
         0.5*displacement(:,2).^2];
    distance = sqrt(sum(displacement.^2,2));
    positiveDistance = distance(distance > 0);
    if isempty(positiveDistance)
        lengthScale = 1;
    else
        lengthScale = median(positiveDistance);
    end
    weight = 1./(1+(distance/max(lengthScale,eps)).^2);
    Aw = bsxfun(@times, A, sqrt(weight));
    bw = sqrt(weight).*psi(patch);
    if rank(Aw) == 6
        coefficient = Aw\bw;
        gradient = coefficient(2:3);
        hessian = [coefficient(4), coefficient(5); ...
                   coefficient(5), coefficient(6)];
        fitResidual = norm(A*coefficient-psi(patch),2) ...
                    / max(norm(psi(patch)-mean(psi(patch)),2),eps);
        eigenvalues = eig(hessian);
        definite = (strcmp(mode,'min') && all(eigenvalues > 0)) || ...
                   (strcmp(mode,'max') && all(eigenvalues < 0));
        if definite && rcond(hessian) > 1e-12
            offset = -hessian\gradient;
            candidate = nodePosition+offset.';
            inside = point_inside_mesh(candidate, nodes, elements);
            patchRadius = max(distance);
            nearby = norm(offset) <= max(1.25*patchRadius, eps);
            if inside && nearby
                fitPosition = candidate;
                fitValue = coefficient(1) + gradient.'*offset ...
                         + 0.5*offset.'*hessian*offset;
                fitAccepted = true;
            end
        end
    end
end

if fitAccepted
    position = fitPosition;
    value = fitValue;
    method = 'local quadratic fit';
else
    position = nodePosition;
    value = nodeValue;
    method = 'interior nodal extremum';
end

axisData = struct();
axisData.R = position(1);
axisData.Z = position(2);
axisData.position = position;
axisData.psi = value;
axisData.mode = mode;
axisData.classification = ['O-' mode];
axisData.method = method;
axisData.refinementAccepted = fitAccepted;
axisData.node = node;
axisData.nodePosition = nodePosition;
axisData.nodePsi = nodeValue;
axisData.boundaryValue = opts.boundaryValue;
axisData.gradientAtNodeFit = gradient;
axisData.hessian = hessian;
axisData.fitResidual = fitResidual;
axisData.patchNodes = patch;

end


function [nodes, elements, boundaryNodes, interiorNodes] = ...
        validate_inputs(meshData, psi)
required = {'nodes','elements','boundaryNodes','interiorNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData, required))
    error('GS:post:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
nodes = meshData.nodes;
elements = meshData.elements;
boundaryNodes = meshData.boundaryNodes(:);
interiorNodes = meshData.interiorNodes(:);
validateattributes(psi, {'numeric'}, ...
    {'real','finite','vector','numel',size(nodes,1)}, mfilename, 'psi');
if isempty(interiorNodes)
    error('GS:post:NoInteriorNodes', 'The mesh has no interior nodes.');
end
end


function opts = validate_options(opts, boundaryValues)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidAxisOptions', 'opts must be a scalar struct.');
end
defaults = struct('mode','auto', ...
    'boundaryValue',median(boundaryValues), 'fitRings',2);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownAxisOption', ...
        'Unknown axis option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.boundaryValue,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.boundaryValue');
validateattributes(opts.fitRings,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1,'<=',6}, ...
    mfilename,'opts.fitRings');
opts.mode = validate_mode(opts.mode);
end


function mode = validate_mode(value)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:post:InvalidAxisMode','mode must be auto, min, or max.');
end
mode = lower(strtrim(value));
if ~ismember(mode,{'auto','min','max'})
    error('GS:post:InvalidAxisMode','mode must be auto, min, or max.');
end
end


function mode = choose_mode(requested, minimum, maximum, boundaryValue)
if ~strcmp(requested,'auto')
    mode = requested;
elseif abs(minimum-boundaryValue) >= abs(maximum-boundaryValue)
    mode = 'min';
else
    mode = 'max';
end
end


function neighbors = vertex_neighbors(elements,nNodes)
edges = [elements(:,[1,2]); elements(:,[2,3]); elements(:,[3,1])];
edges = unique(sort(edges,2),'rows');
neighbors = cell(nNodes,1);
for k = 1:size(edges,1)
    i = edges(k,1);
    j = edges(k,2);
    neighbors{i}(end+1) = j;
    neighbors{j}(end+1) = i;
end
end


function patch = graph_patch(seed,neighbors,nRings)
visited = false(numel(neighbors),1);
visited(seed) = true;
frontier = seed;
for ring = 1:nRings
    next = [];
    for k = 1:numel(frontier)
        next = [next,neighbors{frontier(k)}]; %#ok<AGROW>
    end
    next = unique(next);
    next = next(~visited(next));
    visited(next) = true;
    frontier = next;
    if isempty(frontier)
        break
    end
end
patch = find(visited);
end


function inside = point_inside_mesh(point,nodes,elements)
tri = triangulation(elements,nodes);
inside = ~isnan(pointLocation(tri,point));
end
