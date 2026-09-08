function localEq = build_local_pest_equilibrium(equilibrium,bundle,opts)
%BUILD_LOCAL_PEST_EQUILIBRIUM Construct q and PEST coordinates from FEM GS.
%
%   LOCALEQ = BUILD_LOCAL_PEST_EQUILIBRIUM(EQUILIBRIUM,BUNDLE) evaluates
%   recovered FEM grad(psi), reconstructs F(psi) from the FF' profile, and
%   computes on every surface
%
%       q = (1/2*pi) integral F/(R |grad psi|) dl,
%       V' = 2*pi integral R/|grad psi| dl,
%       sqrt(g)_PEST = q R^2/F.
%
%   The periodic geometric coordinate in BUNDLE is then remapped to the
%   straight-field-line PEST angle. No GEQDSK data are used.
%
%   This function calls, but does not modify, generic functions in the
%   sibling equilibrium folder: compute_contour_integrals,
%   straight_field_line_angles, validate_field_line_straightness, and
%   build_SFL_maps.
%
%   OPTS fields:
%       nTheta       uniform PEST grid size (default bundle.nPoloidal)
%       FBoundary    F=R*Bphi at LCFS (equilibrium cfg value by default)

%   The output contains raw surface data, uniform localEq.map, an FEM field
%   evaluator, q/qprime/sHat/alpha profiles, and numerical diagnostics.

narginchk(2,3);
if nargin<3 || isempty(opts)
    opts = struct();
end
validate_equilibrium(equilibrium);
validate_bundle(bundle);
opts = validate_options(opts,equilibrium,bundle);
add_external_equilibrium_functions();

requiredExternal = {'compute_contour_integrals', ...
    'straight_field_line_angles','validate_field_line_straightness', ...
    'build_SFL_maps','reconstruct_F'};
for k = 1:numel(requiredExternal)
    if exist(requiredExternal{k},'file')~=2
        error('GS:ballooning:MissingEquilibriumFunction', ...
            'Required external function %s is unavailable.', ...
            requiredExternal{k});
    end
end

mesh = equilibrium.mesh;
psi = equilibrium.psi(:);
profiles = equilibrium.profiles;
psiAxis = bundle.psiAxis;
psiBoundary = bundle.psiBoundary;
fluxSpan = psiBoundary-psiAxis;
fluxSpec = struct('psiAxis',psiAxis,'psiBoundary',psiBoundary, ...
    'rangePolicy','clip','tolerance',1e-8);

tri = triangulation(mesh.elements,mesh.nodes);
nodalPsiR = equilibrium.field.node.dpsiDR(:);
nodalPsiZ = equilibrium.field.node.dpsiDZ(:);
FFspec = profiles.FFprime;

eqfunc = struct();
eqfunc.eval = @(R,Z) evaluate_fem_field( ...
    R,Z,tri,mesh.elements,psi,nodalPsiR,nodalPsiZ, ...
    profiles,fluxSpec,opts.FBoundary,FFspec);
eqfunc.psiR = @(R,Z) interpolate_on_mesh( ...
    tri,mesh.elements,nodalPsiR,R,Z);
eqfunc.psiZ = @(R,Z) interpolate_on_mesh( ...
    tri,mesh.elements,nodalPsiZ,R,Z);
eqfunc.F = @(psiN) reconstruct_F( ...
    psiAxis+psiN*fluxSpan,opts.FBoundary,FFspec,fluxSpec);
eqfunc.pprime = @(psiN) pprime( ...
    psiAxis+psiN*fluxSpan,profiles.pprime,fluxSpec);

nSurface = bundle.nSurface;
qIntegral = zeros(nSurface,1);
V = zeros(nSurface,1);
Vprime = zeros(nSurface,1);
pprimeProfile = zeros(nSurface,1);
FProfile = zeros(nSurface,1);
preliminary = repmat(struct(),nSurface,1);
for k = 1:nSurface
    surface = bundle.surfaces(k);
    values = eqfunc.eval(surface.R,surface.Z);
    FProfile(k) = eqfunc.F(surface.psiN);
    pprimeProfile(k) = eqfunc.pprime(surface.psiN);
    qIntegral(k) = sum(0.5*( ...
        FProfile(k)./(surface.R.*values.gradPsi) ...
        +circshift(FProfile(k)./(surface.R.*values.gradPsi),-1)) ...
        .*surface.dl)/(2*pi);
    Vprime(k) = 2*pi*sum(0.5*( ...
        surface.R./values.gradPsi ...
        +circshift(surface.R./values.gradPsi,-1)).*surface.dl);
    [area,intR] = polygon_area_and_Rmoment(surface.R,surface.Z);
    V(k) = 2*pi*abs(intR);
    preliminary(k).field = values;
    preliminary(k).area = abs(area);
end

qprime = local_profile_derivative(bundle.psi,qIntegral);
eqfunc.q = @(psiN) interp1(bundle.psiN,qIntegral,psiN,'pchip','extrap');
eqfunc.qprime = @(psiN) interp1( ...
    bundle.psiN,qprime,psiN,'pchip','extrap');

surfaceData = repmat(struct('surf',[],'ints',[],'ang',[], ...
    'straightness',[],'field',[]),nSurface,1);
for k = 1:nSurface
    base = bundle.surfaces(k);
    surf = struct();
    surf.psiN = base.psiN;
    surf.psi = base.psi;
    surf.R = base.R;
    surf.Z = base.Z;
    surf.area = base.s;
    surf.dl = base.dl;
    surf.L = base.perimeter;
    surf.F = FProfile(k);
    surf.q = qIntegral(k);
    surf.pprime = pprimeProfile(k);
    surf.qprime = qprime(k);
    surf.Raxis = bundle.axis.R;
    surf.Zaxis = bundle.axis.Z;

    % The external function recomputes q, V, V' from the FEM adapter and
    % now uses the integral q profile as its reference rather than GEQDSK.
    ints = compute_contour_integrals(struct(),eqfunc,surf);
    % UseQ="geqdsk" here means the q supplied by eqfunc.q. In this adapter
    % it is not GEQDSK data: it is the segment-trapezoid contour integral
    % computed above, matching the quadrature used by the external PEST
    % angle routine.
    ang = straight_field_line_angles( ...
        eqfunc,surf,ints,UseQ="geqdsk");
    straightness = validate_field_line_straightness(eqfunc,surf,ang);

    values = preliminary(k).field;
    values.JacobianPEST = ints.q_int*surf.R.^2/surf.F;
    values.BthetaContravariant = 1./values.JacobianPEST;
    values.BphiContravariant = ints.q_int./values.JacobianPEST;

    surfaceData(k).surf = surf;
    surfaceData(k).ints = ints;
    surfaceData(k).ang = ang;
    surfaceData(k).straightness = straightness;
    surfaceData(k).field = values;
end

maps = build_SFL_maps(surfaceData,opts.nTheta);
map = maps.PEST;
map.q_profile = qIntegral;
map.qprime = qprime;
map.sHat = 2*V.*qprime./(qIntegral.*Vprime);
map.alpha = arrayfun(@(item) item.ints.alpha,surfaceData);
map.alpha = map.alpha(:);

mapValues = eqfunc.eval(map.R,map.Z);
map.psiEvaluated = mapValues.psi;
map.psiR = mapValues.psiR;
map.psiZ = mapValues.psiZ;
map.gradPsi = mapValues.gradPsi;
map.B2 = mapValues.B2;
map.Bp = mapValues.Bp;
map.Bphi = mapValues.Bphi;
map.Jacobian = repmat(qIntegral,1,opts.nTheta) ...
    .*map.R.^2./repmat(FProfile,1,opts.nTheta);

expectedPsi = repmat(bundle.psi,1,opts.nTheta);
normalizedMapFluxError = (map.psiEvaluated-expectedPsi)/fluxSpan;
straightnessError = arrayfun( ...
    @(item) item.straightness.PEST_linefit_relerr,surfaceData).';
monotone = arrayfun( ...
    @(item) item.ang.diagnostics.PEST_monotone,surfaceData).';
qRelativeError = arrayfun(@(item) item.ints.q_relerr,surfaceData).';

% B^phi=q/J must equal the physical toroidal contravariant F/R^2.
contravariantIdentityError = max(abs( ...
    repmat(qIntegral,1,opts.nTheta)./map.Jacobian ...
    -repmat(FProfile,1,opts.nTheta)./map.R.^2),[],'all');

profilesOut = table(bundle.psiN,bundle.psi,FProfile,qIntegral,qprime, ...
    V,Vprime,map.sHat,map.alpha, ...
    'VariableNames',{'psiN','psi','F','q','qprime','V','Vprime', ...
                     'sHat','alpha'});

localEq = struct();
localEq.bundle = bundle;
localEq.surfaces = surfaceData;
localEq.map = map;
localEq.maps = maps;
localEq.profiles = profilesOut;
localEq.eqfunc = eqfunc;
localEq.fluxSpec = fluxSpec;
localEq.options = opts;
localEq.diagnostics = struct( ...
    'allPESTAnglesMonotone',all(monotone), ...
    'maximumStraightFieldLineError',max(straightnessError), ...
    'maximumQIntegralRelativeError',max(abs(qRelativeError)), ...
    'maximumNormalizedMapFluxError', ...
        max(abs(normalizedMapFluxError),[],'all'), ...
    'rmsNormalizedMapFluxError', ...
        sqrt(mean(normalizedMapFluxError.^2,'all')), ...
    'minimumJacobian',min(map.Jacobian,[],'all'), ...
    'maximumJacobian',max(map.Jacobian,[],'all'), ...
    'contravariantIdentityError',contravariantIdentityError, ...
    'externalFunctionsCalled',{requiredExternal});
end


function values = evaluate_fem_field( ...
        R,Z,tri,elements,psi,psiRNode,psiZNode, ...
        profiles,fluxSpec,FBoundary,FFspec)
psiValue = interpolate_on_mesh(tri,elements,psi,R,Z);
psiR = interpolate_on_mesh(tri,elements,psiRNode,R,Z);
psiZ = interpolate_on_mesh(tri,elements,psiZNode,R,Z);
psiN = normalize_flux(psiValue,fluxSpec.psiAxis, ...
    fluxSpec.psiBoundary,struct('rangePolicy','clip','tolerance',1e-8));
F = reconstruct_F(psiValue,FBoundary,FFspec,fluxSpec);
gradPsi = hypot(psiR,psiZ);
Bp = gradPsi./R;
Bphi = F./R;
B2 = Bp.^2+Bphi.^2;
values = struct('psi',psiValue,'psiN',psiN,'psiR',psiR,'psiZ',psiZ, ...
    'gradPsi',gradPsi,'F',F,'Bp',Bp,'Bphi',Bphi,'B2',B2, ...
    'pprime',pprime(psiValue,profiles.pprime,fluxSpec));
end


function values = interpolate_on_mesh(tri,elements,nodal,R,Z)
shape = size(R);
if ~isequal(size(R),size(Z))
    error('GS:ballooning:CoordinateSizeMismatch', ...
        'R and Z must have equal sizes.');
end
points = [R(:),Z(:)];
triangleIndex = pointLocation(tri,points);
if any(isnan(triangleIndex))
    error('GS:ballooning:PointOutsideFEMMesh', ...
        'A local-equilibrium evaluation point lies outside the FEM mesh.');
end
barycentric = cartesianToBarycentric(tri,triangleIndex,points);
local = reshape(nodal(elements(triangleIndex,:)),size(barycentric));
values = reshape(sum(barycentric.*local,2),shape);
end


function derivative = local_profile_derivative(x,y)
x = x(:);
y = y(:);
n = numel(x);
if n<5
    error('GS:ballooning:InsufficientRadialSurfaces', ...
        'At least five surfaces are required for local radial derivatives.');
end
derivative = zeros(n,1);
for k = 1:n
    [~,order] = sort(abs(x-x(k)));
    index = sort(order(1:5));
    dx = x(index)-x(k);
    degree = min(3,numel(index)-1);
    matrix = zeros(numel(index),degree+1);
    for power = 0:degree
        matrix(:,power+1) = dx.^power;
    end
    coefficient = matrix\y(index);
    derivative(k) = coefficient(2);
end
end


function [area,intR] = polygon_area_and_Rmoment(R,Z)
nextR = R([2:end,1]);
nextZ = Z([2:end,1]);
cross = R.*nextZ-nextR.*Z;
area = 0.5*sum(cross);
intR = (1/6)*sum((R+nextR).*cross);
end


function validate_equilibrium(equilibrium)
required = {'mesh','psi','axis','profiles','field','configuration'};
if ~isstruct(equilibrium) || ~isscalar(equilibrium) || ...
        ~all(isfield(equilibrium,required)) || ...
        ~isfield(equilibrium.field,'node') || ...
        ~all(isfield(equilibrium.field.node,{'dpsiDR','dpsiDZ'}))
    error('GS:ballooning:InvalidEquilibrium', ...
        'equilibrium must be the completed result from main_cos.');
end
end


function validate_bundle(bundle)
required = {'surfaces','psi','psiN','psiAxis','psiBoundary', ...
    'nSurface','nPoloidal','axis'};
if ~isstruct(bundle) || ~isscalar(bundle) || ...
        ~all(isfield(bundle,required)) || bundle.nSurface<5
    error('GS:ballooning:InvalidSurfaceBundle', ...
        'Use a bundle of at least five surfaces from stage 1.');
end
end


function opts = validate_options(opts,equilibrium,bundle)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:ballooning:InvalidLocalEquilibriumOptions', ...
        'opts must be a scalar struct.');
end
if isfield(equilibrium.configuration,'FBoundary')
    defaultF = equilibrium.configuration.FBoundary;
else
    defaultF = 3;
end
defaults = struct('nTheta',bundle.nPoloidal,'FBoundary',defaultF);
unknown = setdiff(fieldnames(opts),fieldnames(defaults));
if ~isempty(unknown)
    error('GS:ballooning:UnknownLocalEquilibriumOption', ...
        'Unknown option: %s',strjoin(unknown,', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts,names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
validateattributes(opts.nTheta,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',64},mfilename,'opts.nTheta');
validateattributes(opts.FBoundary,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'opts.FBoundary');
end


function add_external_equilibrium_functions()
gsRoot = fileparts(fileparts(mfilename('fullpath')));
equilibriumFolder = fullfile(fileparts(gsRoot),'equilibrium');
if ~isfolder(equilibriumFolder)
    error('GS:ballooning:MissingEquilibriumFolder', ...
        'The sibling equilibrium function folder was not found.');
end
addpath(equilibriumFolder);
end
