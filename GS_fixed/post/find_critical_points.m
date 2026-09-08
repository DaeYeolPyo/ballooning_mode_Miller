function critical = find_critical_points(meshData,psi,opts)
%FIND_CRITICAL_POINTS Find discrete interior O- and X-point candidates.
%
%   CRITICAL = FIND_CRITICAL_POINTS(MESH,PSI) analyzes the cyclic signs of
%   psi(neighbor)-psi(node) around every interior mesh vertex. Local extrema
%   are O-point candidates; four or more sign transitions are PL saddle
%   (X-point) candidates.
%
%   This is a topology diagnostic for a P1 field, not a sub-element Newton
%   search. The primary O-point should be refined separately with find_axis.

narginchk(2,3);
required = {'nodes','elements','interiorNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData,required))
    error('GS:post:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
nodes = meshData.nodes;
elements = meshData.elements;
interiorNodes = meshData.interiorNodes(:);
validateattributes(psi,{'numeric'}, ...
    {'real','finite','vector','numel',size(nodes,1)},mfilename,'psi');
psi = psi(:);

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts,psi);
neighbors = ordered_neighbors(nodes,elements);

critical = struct('node',{},'R',{},'Z',{},'psi',{}, ...
    'type',{},'signTransitions',{},'prominence',{});
for k = 1:numel(interiorNodes)
    node = interiorNodes(k);
    ring = neighbors{node};
    if numel(ring) < 3
        continue
    end
    difference = psi(ring)-psi(node);
    prominence = min(max(difference),-min(difference));
    if all(difference >= -opts.valueTolerance) && ...
            any(difference > opts.valueTolerance)
        type = 'O-min';
        transitions = 0;
        prominence = max(difference);
    elseif all(difference <= opts.valueTolerance) && ...
            any(difference < -opts.valueTolerance)
        type = 'O-max';
        transitions = 0;
        prominence = -min(difference);
    else
        signs = sign(difference);
        signs(abs(difference)<=opts.valueTolerance) = 0;
        signs = signs(signs~=0);
        if numel(signs) < 4
            continue
        end
        transitions = nnz(signs~=circshift(signs,1));
        if transitions < 4
            continue
        end
        type = 'X';
    end
    if prominence < opts.minimumProminence
        continue
    end
    point = struct('node',node,'R',nodes(node,1),'Z',nodes(node,2), ...
        'psi',psi(node),'type',type,'signTransitions',transitions, ...
        'prominence',prominence);
    critical(end+1,1) = point; %#ok<AGROW>
end

end


function neighbors = ordered_neighbors(nodes,elements)
nNodes = size(nodes,1);
edges = unique(sort([elements(:,[1,2]);elements(:,[2,3]); ...
    elements(:,[3,1])],2),'rows');
neighbors = cell(nNodes,1);
for k = 1:size(edges,1)
    i = edges(k,1);
    j = edges(k,2);
    neighbors{i}(end+1) = j;
    neighbors{j}(end+1) = i;
end
for node = 1:nNodes
    ring = neighbors{node};
    if isempty(ring)
        continue
    end
    delta = nodes(ring,:)-nodes(node,:);
    [~,order] = sort(atan2(delta(:,2),delta(:,1)));
    neighbors{node} = ring(order);
end
end


function opts = validate_options(opts,psi)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidCriticalOptions','opts must be a scalar struct.');
end
scale = max(1,max(psi)-min(psi));
defaults = struct('valueTolerance',1e-10*scale, ...
                  'minimumProminence',1e-8*scale);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownCriticalOption', ...
        'Unknown critical-point option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.valueTolerance,{'numeric'}, ...
    {'real','finite','scalar','nonnegative'},mfilename,'opts.valueTolerance');
validateattributes(opts.minimumProminence,{'numeric'}, ...
    {'real','finite','scalar','nonnegative'}, ...
    mfilename,'opts.minimumProminence');
end
