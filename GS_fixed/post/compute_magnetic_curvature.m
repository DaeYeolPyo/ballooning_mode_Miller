function result = compute_magnetic_curvature( ...
    meshData,psi,field,profiles,opts)
%COMPUTE_MAGNETIC_CURVATURE Magnetic curvature and kappa dot grad(p).
%
%   RESULT = COMPUTE_MAGNETIC_CURVATURE(MESH,PSI,FIELD,PROFILES,OPTS)
%   evaluates the axisymmetric cylindrical-coordinate curvature
%
%       kappa = (b dot grad)b,          b = B/|B|,
%
%   and the signed pressure-curvature measure
%
%       kappa dot grad(p) = p'(psi) kappa dot grad(psi).
%
%   FIELD must be returned by compute_B with a nonempty toroidal field.
%   PROFILES is returned by make_profiles. OPTS supports:
%       psiAxis             flux at the magnetic axis (inferred if empty)
%       boundaryValue       LCFS flux (boundary median by default)
%       projectPerpendicular remove numerical b-parallel curvature (true)
%       normalizedFluxEdges edges for volume-weighted band statistics
%
%   Axisymmetry does not eliminate cylindrical basis-vector derivatives.
%   The element formulas are
%
%       kappa_R   = b_R*d_R b_R + b_Z*d_Z b_R - b_phi^2/R,
%       kappa_phi = b_R*d_R b_phi + b_Z*d_Z b_phi + b_R*b_phi/R,
%       kappa_Z   = b_R*d_R b_Z + b_Z*d_Z b_Z.
%
%   Derivatives use a recovered nodal b followed by a P1 element gradient.
%   Element values are primary; nodal values are area-weighted recoveries.

narginchk(4,5);
if nargin < 5 || isempty(opts)
    opts = struct();
end
[nodes,elements,area,boundaryNodes] = validate_mesh(meshData);
nNodes = size(nodes,1);
nElements = size(elements,1);
validateattributes(psi,{'numeric'}, ...
    {'real','finite','vector','numel',nNodes},mfilename,'psi');
psi = psi(:);
validate_field(field,nNodes,nElements);
validate_profiles(profiles);
opts = validate_options(opts,psi,boundaryNodes);

% Dynamic solver axis values are not written back into cfg.profiles, so the
% diagnostic explicitly constructs the final-equilibrium flux convention.
fluxSpec = struct('psiAxis',opts.psiAxis, ...
    'psiBoundary',opts.boundaryValue, ...
    'rangePolicy','clip','tolerance',1e-8);

nodeBtotal = field.node.Btotal(:);
elementBtotal = field.element.Btotal(:);
fieldScale = max([1; abs(nodeBtotal); abs(elementBtotal)]);
fieldTol = 100*eps(max(fieldScale));
if any(nodeBtotal<=fieldTol) || any(elementBtotal<=fieldTol)
    error('GS:post:VanishingMagneticField', ...
        ['Magnetic curvature is undefined where |B| vanishes. Supply a ' ...
         'nonzero toroidal field to compute_B and inspect the equilibrium.']);
end

% Recovered nodal unit vector is differentiated within each P1 element.
nodeBR = field.node.BR(:)./nodeBtotal;
nodeBphi = field.node.Bphi(:)./nodeBtotal;
nodeBZ = field.node.BZ(:)./nodeBtotal;
[dbRDR,dbRDZ] = p1_gradient(nodes,elements,nodeBR);
[dbPhiDR,dbPhiDZ] = p1_gradient(nodes,elements,nodeBphi);
[dbZDR,dbZDZ] = p1_gradient(nodes,elements,nodeBZ);

elementBR = field.element.BR(:)./elementBtotal;
elementBphi = field.element.Bphi(:)./elementBtotal;
elementBZ = field.element.BZ(:)./elementBtotal;
centroid = field.element.centroid;
R = centroid(:,1);

kappaRRaw = elementBR.*dbRDR + elementBZ.*dbRDZ ...
          - elementBphi.^2./R;
kappaPhiRaw = elementBR.*dbPhiDR + elementBZ.*dbPhiDZ ...
            + elementBR.*elementBphi./R;
kappaZRaw = elementBR.*dbZDR + elementBZ.*dbZDZ;
parallelRaw = kappaRRaw.*elementBR ...
            + kappaPhiRaw.*elementBphi + kappaZRaw.*elementBZ;

if opts.projectPerpendicular
    kappaR = kappaRRaw-parallelRaw.*elementBR;
    kappaPhi = kappaPhiRaw-parallelRaw.*elementBphi;
    kappaZ = kappaZRaw-parallelRaw.*elementBZ;
else
    kappaR = kappaRRaw;
    kappaPhi = kappaPhiRaw;
    kappaZ = kappaZRaw;
end
kappaMagnitude = sqrt(kappaR.^2+kappaPhi.^2+kappaZ.^2);
parallelAfterProjection = kappaR.*elementBR ...
    +kappaPhi.*elementBphi+kappaZ.*elementBZ;

% Pressure is a flux function, so grad(p)=p'(psi)*grad(psi).
[dpsiDR,dpsiDZ] = p1_gradient(nodes,elements,psi);
elementPsi = mean(reshape(psi(elements(:)),size(elements)),2);
elementPprime = pprime(elementPsi,profiles.pprime,fluxSpec);
gradPR = elementPprime.*dpsiDR;
gradPZ = elementPprime.*dpsiDZ;
kappaDotGradP = kappaR.*gradPR+kappaZ.*gradPZ;

% Nodal visualization fields use conservative area-weighted recovery.
nodeKappaR = recover_to_nodes(elements,area,kappaR,nNodes);
nodeKappaPhi = recover_to_nodes(elements,area,kappaPhi,nNodes);
nodeKappaZ = recover_to_nodes(elements,area,kappaZ,nNodes);
nodeKappaMagnitude = sqrt( ...
    nodeKappaR.^2+nodeKappaPhi.^2+nodeKappaZ.^2);
nodeGradPsiR = recover_to_nodes(elements,area,dpsiDR,nNodes);
nodeGradPsiZ = recover_to_nodes(elements,area,dpsiDZ,nNodes);
nodePprime = pprime(psi,profiles.pprime,fluxSpec);
nodeGradPR = nodePprime.*nodeGradPsiR;
nodeGradPZ = nodePprime.*nodeGradPsiZ;
nodeKappaDotGradP = nodeKappaR.*nodeGradPR+nodeKappaZ.*nodeGradPZ;

elementPsiNormalized = normalize_flux( ...
    elementPsi,opts.psiAxis,opts.boundaryValue, ...
    struct('rangePolicy','clip','tolerance',1e-8));
nodePsiNormalized = normalize_flux( ...
    psi,opts.psiAxis,opts.boundaryValue, ...
    struct('rangePolicy','clip','tolerance',1e-8));

element = struct();
element.centroid = centroid;
element.psi = elementPsi;
element.psiNormalized = elementPsiNormalized;
element.bR = elementBR;
element.bphi = elementBphi;
element.bZ = elementBZ;
element.kappaRRaw = kappaRRaw;
element.kappaPhiRaw = kappaPhiRaw;
element.kappaZRaw = kappaZRaw;
element.kappaParallelRaw = parallelRaw;
element.kappaR = kappaR;
element.kappaPhi = kappaPhi;
element.kappaZ = kappaZ;
element.kappaMagnitude = kappaMagnitude;
element.kappaParallel = parallelAfterProjection;
element.pprime = elementPprime;
element.gradPsiR = dpsiDR;
element.gradPsiZ = dpsiDZ;
element.gradPR = gradPR;
element.gradPZ = gradPZ;
element.kappaDotGradP = kappaDotGradP;

node = struct();
node.position = nodes;
node.psi = psi;
node.psiNormalized = nodePsiNormalized;
node.bR = nodeBR;
node.bphi = nodeBphi;
node.bZ = nodeBZ;
node.kappaR = nodeKappaR;
node.kappaPhi = nodeKappaPhi;
node.kappaZ = nodeKappaZ;
node.kappaMagnitude = nodeKappaMagnitude;
node.pprime = nodePprime;
node.gradPsiR = nodeGradPsiR;
node.gradPsiZ = nodeGradPsiZ;
node.gradPR = nodeGradPR;
node.gradPZ = nodeGradPZ;
node.kappaDotGradP = nodeKappaDotGradP;

stats = curvature_statistics( ...
    R,area,elementPsiNormalized,kappaMagnitude,kappaDotGradP, ...
    opts.normalizedFluxEdges);
stats.maximumRawParallelFraction = max(abs(parallelRaw) ...
    ./max(sqrt(kappaRRaw.^2+kappaPhiRaw.^2+kappaZRaw.^2),eps));
stats.maximumProjectedParallelFraction = max(abs(parallelAfterProjection) ...
    ./max(kappaMagnitude,eps));

result = struct();
result.element = element;
result.node = node;
result.statistics = stats;
result.flux = fluxSpec;
result.options = opts;
result.convention = struct( ...
    'coordinates','axisymmetric cylindrical (R,phi,Z)', ...
    'pressureGradient','grad(p)=pprime(psi)*grad(psi)', ...
    'signedMeasure','kappa dot grad(p)', ...
    'signNote',[ ...
        'The result is signed and no sign is relabeled automatically. ' ...
        'Interpret good/bad curvature using the chosen pressure and flux ' ...
        'orientation.']);

end


function stats = curvature_statistics( ...
        R,area,psiN,kappaMagnitude,measure,edges)
areaWeight = area;
volumeWeight = 2*pi*R.*area;
measureScale = max([1;abs(measure)]);
zeroTolerance = 1e-12*max(measureScale);
positive = measure>zeroTolerance;
negative = measure<-zeroTolerance;
nearZero = ~(positive|negative);

stats = struct();
stats.minimumKappaMagnitude = min(kappaMagnitude);
stats.maximumKappaMagnitude = max(kappaMagnitude);
stats.areaWeightedMeanKappaMagnitude = ...
    sum(areaWeight.*kappaMagnitude)/sum(areaWeight);
stats.minimumKappaDotGradP = min(measure);
stats.maximumKappaDotGradP = max(measure);
stats.areaWeightedMeanKappaDotGradP = ...
    sum(areaWeight.*measure)/sum(areaWeight);
stats.volumeWeightedMeanKappaDotGradP = ...
    sum(volumeWeight.*measure)/sum(volumeWeight);
stats.positiveAreaFraction = sum(areaWeight(positive))/sum(areaWeight);
stats.negativeAreaFraction = sum(areaWeight(negative))/sum(areaWeight);
stats.nearZeroAreaFraction = sum(areaWeight(nearZero))/sum(areaWeight);
stats.positiveVolumeFraction = ...
    sum(volumeWeight(positive))/sum(volumeWeight);
stats.negativeVolumeFraction = ...
    sum(volumeWeight(negative))/sum(volumeWeight);
stats.nearZeroVolumeFraction = ...
    sum(volumeWeight(nearZero))/sum(volumeWeight);
activeVolumeFraction = stats.positiveVolumeFraction ...
    +stats.negativeVolumeFraction;
stats.positiveActiveVolumeFraction = ...
    stats.positiveVolumeFraction/activeVolumeFraction;
stats.negativeActiveVolumeFraction = ...
    stats.negativeVolumeFraction/activeVolumeFraction;
stats.signTolerance = zeroTolerance;

nBands = numel(edges)-1;
lower = edges(1:end-1);
upper = edges(2:end);
nElements = zeros(nBands,1);
meanValue = NaN(nBands,1);
minimumValue = NaN(nBands,1);
maximumValue = NaN(nBands,1);
positiveVolumeFraction = NaN(nBands,1);
for k = 1:nBands
    if k<nBands
        inBand = psiN>=lower(k) & psiN<upper(k);
    else
        inBand = psiN>=lower(k) & psiN<=upper(k);
    end
    nElements(k) = nnz(inBand);
    if ~any(inBand)
        continue
    end
    weights = volumeWeight(inBand);
    values = measure(inBand);
    meanValue(k) = sum(weights.*values)/sum(weights);
    minimumValue(k) = min(values);
    maximumValue(k) = max(values);
    positiveVolumeFraction(k) = ...
        sum(weights(values>zeroTolerance))/sum(weights);
end
stats.fluxBands = table(lower,upper,nElements,meanValue, ...
    minimumValue,maximumValue,positiveVolumeFraction);
end


function nodal = recover_to_nodes(elements,area,values,nNodes)
weights = accumarray(elements(:),repmat(area,3,1), ...
    [nNodes,1],@sum,0);
nodal = accumarray(elements(:),repmat(area.*values,3,1), ...
    [nNodes,1],@sum,0)./weights;
end


function [gradientR,gradientZ] = p1_gradient(nodes,elements,values)
p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
twiceArea = (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) ...
    -(p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1));
gradNR = [(p2(:,2)-p3(:,2))./twiceArea, ...
          (p3(:,2)-p1(:,2))./twiceArea, ...
          (p1(:,2)-p2(:,2))./twiceArea];
gradNZ = [(p3(:,1)-p2(:,1))./twiceArea, ...
          (p1(:,1)-p3(:,1))./twiceArea, ...
          (p2(:,1)-p1(:,1))./twiceArea];
elementValues = reshape(values(elements(:)),size(elements));
gradientR = sum(elementValues.*gradNR,2);
gradientZ = sum(elementValues.*gradNZ,2);
end


function [nodes,elements,area,boundaryNodes] = validate_mesh(meshData)
required = {'nodes','elements','elementArea','boundaryNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData,required))
    error('GS:post:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
nodes = meshData.nodes;
elements = meshData.elements;
area = meshData.elementArea(:);
boundaryNodes = meshData.boundaryNodes(:);
if size(nodes,2)~=2 || any(~isfinite(nodes(:))) || any(nodes(:,1)<=0) || ...
        size(elements,2)~=3 || numel(area)~=size(elements,1) || any(area<=0)
    error('GS:post:InvalidMesh','The mesh arrays are inconsistent.');
end
end


function validate_field(field,nNodes,nElements)
locations = {'node','element'};
counts = [nNodes,nElements];
for k = 1:2
    location = locations{k};
    if ~isfield(field,location) || ~isstruct(field.(location))
        error('GS:post:InvalidMagneticField', ...
            'field must be returned by compute_B.');
    end
    required = {'BR','Bphi','BZ','Btotal'};
    data = field.(location);
    if ~all(isfield(data,required))
        error('GS:post:MissingToroidalField', ...
            'compute_B must be called with F so Bphi and Btotal exist.');
    end
    for j = 1:numel(required)
        value = data.(required{j});
        if ~isnumeric(value) || numel(value)~=counts(k) || ...
                ~isreal(value) || any(~isfinite(value(:)))
            error('GS:post:InvalidMagneticField', ...
                'field.%s.%s is invalid.',location,required{j});
        end
    end
end
if ~isfield(field.element,'centroid') || ...
        ~isequal(size(field.element.centroid),[nElements,2])
    error('GS:post:InvalidMagneticField', ...
        'field.element.centroid is missing or inconsistent.');
end
end


function validate_profiles(profiles)
if ~isstruct(profiles) || ~isscalar(profiles) || ...
        ~isfield(profiles,'pprime')
    error('GS:post:InvalidProfiles', ...
        'profiles must contain a pprime specification from make_profiles.');
end
if exist('pprime','file')~=2 || exist('normalize_flux','file')~=2
    error('GS:post:MissingPhysicsDependency', ...
        'Add the physics folder to the MATLAB path.');
end
end


function opts = validate_options(opts,psi,boundaryNodes)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:post:InvalidCurvatureOptions','opts must be a scalar struct.');
end
boundaryDefault = median(psi(boundaryNodes));
distance = abs(psi-boundaryDefault);
[~,axisNode] = max(distance);
defaults = struct('psiAxis',psi(axisNode), ...
    'boundaryValue',boundaryDefault,'projectPerpendicular',true, ...
    'normalizedFluxEdges',(0:0.1:1).');
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:post:UnknownCurvatureOption', ...
        'Unknown curvature option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k}) || isempty(opts.(names{k}))
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.psiAxis,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.psiAxis');
validateattributes(opts.boundaryValue,{'numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.boundaryValue');
if opts.psiAxis==opts.boundaryValue
    error('GS:post:CollapsedFluxRange', ...
        'psiAxis and boundaryValue must differ.');
end
validateattributes(opts.projectPerpendicular,{'logical','numeric'}, ...
    {'real','finite','scalar'},mfilename,'opts.projectPerpendicular');
opts.projectPerpendicular = logical(opts.projectPerpendicular);
validateattributes(opts.normalizedFluxEdges,{'numeric'}, ...
    {'real','finite','vector','>=',0,'<=',1}, ...
    mfilename,'opts.normalizedFluxEdges');
opts.normalizedFluxEdges = opts.normalizedFluxEdges(:);
if numel(opts.normalizedFluxEdges)<2 || ...
        opts.normalizedFluxEdges(1)~=0 || ...
        opts.normalizedFluxEdges(end)~=1 || ...
        any(diff(opts.normalizedFluxEdges)<=0)
    error('GS:post:InvalidFluxEdges', ...
        'normalizedFluxEdges must increase strictly from 0 to 1.');
end
end
