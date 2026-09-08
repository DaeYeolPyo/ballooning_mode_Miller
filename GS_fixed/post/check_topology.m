function topology = check_topology(meshData,psi,axisData,opts)
%CHECK_TOPOLOGY Check nested closed flux contours and interior critical points.
%
%   TOPOLOGY = CHECK_TOPOLOGY(MESH,PSI) locates the axis and checks the
%   normalized levels 0.1:0.1:0.9. A valid single-axis topology requires at
%   each level exactly one closed contour, no open contour, and that contour
%   must enclose the same axis. Consecutive contours are also checked for
%   geometric nesting.
%
%   TOPOLOGY = CHECK_TOPOLOGY(MESH,PSI,AXISDATA,OPTS) accepts:
%       normalizedLevels  levels strictly between 0 and 1
%       boundaryValue     scalar LCFS flux (boundary median by default)
%       failOnFailure     throw when checks fail (default false)
%       requireSingleO    require one discrete O candidate (default true)
%       requireNoX        require no discrete X candidate (default true)

%   Critical-point detection is a P1 vertex diagnostic; contour checks are
%   the primary topology criterion.

narginchk(2,4);
required = {'nodes','boundaryNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData,required))
    error('GS:post:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
validateattributes(psi,{'numeric'}, ...
    {'real','finite','vector','numel',size(meshData.nodes,1)},mfilename,'psi');
psi = psi(:);

if nargin < 3 || isempty(axisData)
    axisData = [];
end
if nargin < 4 || isempty(opts)
    opts = struct();
end
boundaryDefault = median(psi(meshData.boundaryNodes));
opts = validate_options(opts,boundaryDefault);
if isempty(axisData)
    axisData = find_axis(meshData,psi, ...
        struct('boundaryValue',opts.boundaryValue));
end
validate_axis(axisData);

fluxSpan = opts.boundaryValue-axisData.psi;
if abs(fluxSpan)<=100*eps(max([1,abs(opts.boundaryValue),abs(axisData.psi)]))
    error('GS:post:CollapsedFluxRange', ...
        'Axis and boundary flux values are numerically indistinguishable.');
end
levels = axisData.psi+opts.normalizedLevels*fluxSpan;
contourSets = extract_flux_contours(meshData,psi,levels);

nLevels = numel(levels);
levelChecks = repmat(struct('psiNormalized',[],'psi',[], ...
    'nCurves',0,'nClosed',0,'nOpen',0,'nEnclosingAxis',0, ...
    'singleClosedAxisContour',false,'enclosingCurve',[]),nLevels,1);
selectedCurves = cell(nLevels,1);
for k = 1:nLevels
    curves = contourSets(k).curves;
    contains = false(numel(curves),1);
    for j = 1:numel(curves)
        if curves(j).closed
            points = curves(j).points;
            contains(j) = inpolygon(axisData.R,axisData.Z, ...
                points(:,1),points(:,2));
        end
    end
    enclosing = find(contains);
    valid = contourSets(k).nClosed==1 && ...
            contourSets(k).nOpen==0 && numel(enclosing)==1;
    if numel(enclosing)==1
        selectedCurves{k} = curves(enclosing).points;
    else
        selectedCurves{k} = zeros(0,2);
    end
    levelChecks(k).psiNormalized = opts.normalizedLevels(k);
    levelChecks(k).psi = levels(k);
    levelChecks(k).nCurves = contourSets(k).nCurves;
    levelChecks(k).nClosed = contourSets(k).nClosed;
    levelChecks(k).nOpen = contourSets(k).nOpen;
    levelChecks(k).nEnclosingAxis = numel(enclosing);
    levelChecks(k).singleClosedAxisContour = valid;
    levelChecks(k).enclosingCurve = selectedCurves{k};
end

nestedPairs = true(max(nLevels-1,0),1);
for k = 1:nLevels-1
    inner = selectedCurves{k};
    outer = selectedCurves{k+1};
    if isempty(inner) || isempty(outer)
        nestedPairs(k) = false;
        continue
    end
    if size(inner,1)>1 && isequal(inner(1,:),inner(end,:))
        inner = inner(1:end-1,:);
    end
    [inside,on] = inpolygon(inner(:,1),inner(:,2), ...
        outer(:,1),outer(:,2));
    nestedPairs(k) = all(inside|on);
end

critical = find_critical_points(meshData,psi);
types = {critical.type};
nO = nnz(startsWith(types,'O-'));
nX = nnz(strcmp(types,'X'));
contourPassed = all([levelChecks.singleClosedAxisContour]) ...
             && all(nestedPairs);
criticalPassed = (~opts.requireSingleO || nO==1) ...
              && (~opts.requireNoX || nX==0);
passed = contourPassed && criticalPassed;

topology = struct();
topology.passed = passed;
topology.contourPassed = contourPassed;
topology.criticalPassed = criticalPassed;
topology.axis = axisData;
topology.boundaryValue = opts.boundaryValue;
topology.normalizedLevels = opts.normalizedLevels;
topology.levels = levels;
topology.contourSets = contourSets;
topology.levelChecks = levelChecks;
topology.nestedPairs = nestedPairs;
topology.criticalPoints = critical;
topology.nOPoints = nO;
topology.nXPoints = nX;
topology.options = opts;

if ~passed && opts.failOnFailure
    error('GS:post:InvalidFluxTopology', ...
        ['Flux topology check failed: contours=%d, critical points=%d ' ...
         '(O=%d, X=%d).'],contourPassed,criticalPassed,nO,nX);
end

end


function opts = validate_options(opts,boundaryDefault)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidTopologyOptions','opts must be a scalar struct.');
end
defaults = struct('normalizedLevels',(0.1:0.1:0.9).', ...
    'boundaryValue',boundaryDefault,'failOnFailure',false, ...
    'requireSingleO',true,'requireNoX',true);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownTopologyOption', ...
        'Unknown topology option: %s',strjoin(unknown,', '));
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
if any(diff(opts.normalizedLevels)<=0)
    error('GS:post:InvalidTopologyLevels', ...
        'normalizedLevels must be strictly increasing.');
end
validateattributes(opts.boundaryValue,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.boundaryValue');
opts.failOnFailure = logical_scalar(opts.failOnFailure,'failOnFailure');
opts.requireSingleO = logical_scalar(opts.requireSingleO,'requireSingleO');
opts.requireNoX = logical_scalar(opts.requireNoX,'requireNoX');
end


function value = logical_scalar(value,name)
validateattributes(value,{'logical','numeric'}, ...
    {'real','finite','scalar'},mfilename,['opts.' name]);
value = logical(value);
end


function validate_axis(axisData)
required = {'R','Z','psi'};
if ~isstruct(axisData) || ~isscalar(axisData) || ...
        ~all(isfield(axisData,required)) || ...
        any(~isfinite([axisData.R,axisData.Z,axisData.psi]))
    error('GS:post:InvalidAxisData', ...
        'axisData must be a finite scalar result from find_axis.');
end
end
