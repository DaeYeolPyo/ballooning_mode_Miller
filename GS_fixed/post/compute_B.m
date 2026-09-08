function field = compute_B(meshData,psi,F)
%COMPUTE_B Compute axisymmetric magnetic-field components from a P1 psi.
%
%   FIELD = COMPUTE_B(MESH,PSI) computes
%
%       B_R = -(1/R)*dpsi/dZ,   B_Z = (1/R)*dpsi/dR
%
%   exactly within every P1 element and by area-weighted recovery at nodes.
%
%   FIELD = COMPUTE_B(MESH,PSI,F) additionally computes B_phi=F/R. F may
%   be a scalar, one nodal value per mesh node, or a callback with signature
%   F(psi), F(psi,R), or F(psi,R,Z). If F is omitted, Bphi fields are empty.
%
%   FIELD.element and FIELD.node contain positions, flux gradients, BR, BZ,
%   Bp, and optional Bphi/Btotal values.

narginchk(2,3);
[nodes,elements,area] = validate_inputs(meshData,psi);
psi = psi(:);
if nargin < 3
    F = [];
end

p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
twiceArea = 2*area;
elementPsi = reshape(psi(elements(:)),size(elements));

gradNR = [(p2(:,2)-p3(:,2))./twiceArea, ...
          (p3(:,2)-p1(:,2))./twiceArea, ...
          (p1(:,2)-p2(:,2))./twiceArea];
gradNZ = [(p3(:,1)-p2(:,1))./twiceArea, ...
          (p1(:,1)-p3(:,1))./twiceArea, ...
          (p2(:,1)-p1(:,1))./twiceArea];
dpsiDR = sum(elementPsi.*gradNR,2);
dpsiDZ = sum(elementPsi.*gradNZ,2);
centroid = (p1+p2+p3)/3;
elementBR = -dpsiDZ./centroid(:,1);
elementBZ = dpsiDR./centroid(:,1);
elementBp = hypot(elementBR,elementBZ);

nNodes = size(nodes,1);
weightSum = accumarray(elements(:),repmat(area,3,1),[nNodes,1],@sum,0);
nodalDpsiDR = accumarray(elements(:), ...
    repmat(area.*dpsiDR,3,1),[nNodes,1],@sum,0)./weightSum;
nodalDpsiDZ = accumarray(elements(:), ...
    repmat(area.*dpsiDZ,3,1),[nNodes,1],@sum,0)./weightSum;
nodalBR = -nodalDpsiDZ./nodes(:,1);
nodalBZ = nodalDpsiDR./nodes(:,1);
nodalBp = hypot(nodalBR,nodalBZ);

if isempty(F)
    nodalF = [];
    elementF = [];
    nodalBphi = [];
    elementBphi = [];
    nodalBtotal = [];
    elementBtotal = [];
else
    nodalF = evaluate_F(F,psi,nodes(:,1),nodes(:,2),nNodes,'node');
    elementPsiCentroid = mean(elementPsi,2);
    elementF = evaluate_element_F( ...
        F,nodalF,elements,elementPsiCentroid,centroid);
    nodalBphi = nodalF./nodes(:,1);
    elementBphi = elementF./centroid(:,1);
    nodalBtotal = hypot(nodalBp,nodalBphi);
    elementBtotal = hypot(elementBp,elementBphi);
end

field = struct();
field.node = struct( ...
    'position',nodes,'dpsiDR',nodalDpsiDR,'dpsiDZ',nodalDpsiDZ, ...
    'BR',nodalBR,'BZ',nodalBZ,'Bp',nodalBp, ...
    'F',nodalF,'Bphi',nodalBphi,'Btotal',nodalBtotal);
field.element = struct( ...
    'centroid',centroid,'area',area, ...
    'dpsiDR',dpsiDR,'dpsiDZ',dpsiDZ, ...
    'BR',elementBR,'BZ',elementBZ,'Bp',elementBp, ...
    'F',elementF,'Bphi',elementBphi,'Btotal',elementBtotal);

end


function values = evaluate_element_F(F,nodalF,elements,psi,position)
if isnumeric(F) && ~isscalar(F)
    elementNodalF = reshape(nodalF(elements(:)),size(elements));
    values = mean(elementNodalF,2);
else
    values = evaluate_F(F,psi,position(:,1),position(:,2), ...
        size(position,1),'element');
end
end


function values = evaluate_F(F,psi,R,Z,nExpected,location)
if isnumeric(F)
    if isscalar(F)
        values = repmat(F,nExpected,1);
    elseif isvector(F) && numel(F)==nExpected
        values = F(:);
    else
        error('GS:post:InvalidFValues', ...
            'Numeric F must be scalar or one value per %s.',location);
    end
elseif isa(F,'function_handle')
    nInputs = nargin(F);
    if nInputs == 1
        values = F(psi);
    elseif nInputs == 2
        values = F(psi,R);
    else
        values = F(psi,R,Z);
    end
    if isscalar(values)
        values = repmat(values,nExpected,1);
    else
        values = values(:);
    end
else
    error('GS:post:InvalidFSpecification', ...
        'F must be numeric or a function handle.');
end
if numel(values)~=nExpected || ~isreal(values) || ...
        any(~isfinite(values))
    error('GS:post:InvalidFOutput', ...
        'F evaluation must return one finite real value per %s.',location);
end
end


function [nodes,elements,area] = validate_inputs(meshData,psi)
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
if any(nodes(:,1)<=0)
    error('GS:post:NonPositiveRadius','All mesh nodes must have R>0.');
end
p1 = nodes(elements(:,1),:);
p2 = nodes(elements(:,2),:);
p3 = nodes(elements(:,3),:);
twiceArea = (p2(:,1)-p1(:,1)).*(p3(:,2)-p1(:,2)) ...
          -(p2(:,2)-p1(:,2)).*(p3(:,1)-p1(:,1));
if any(twiceArea<=0)
    error('GS:post:InvalidElementOrientation', ...
        'All mesh triangles must be counterclockwise and nondegenerate.');
end
area = 0.5*twiceArea;
end
