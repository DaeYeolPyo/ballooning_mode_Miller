function report = verify_post(doPlot)
%VERIFY_POST Verify magnetic axis, contours, topology, and magnetic field.

if nargin < 1
    doPlot = true;
end
validateattributes(doPlot,{'logical','numeric'}, ...
    {'real','finite','scalar'},mfilename,'doPlot');
doPlot = logical(doPlot);

dependencies = {'make_boundary','generate_mesh','assemble_stiffness', ...
    'initial_guess','make_profiles','solve_picard'};
for k = 1:numel(dependencies)
    if exist(dependencies{k},'file')~=2
        error('GS:post:MissingDependency', ...
            'Required function %s is not on the MATLAB path.',dependencies{k});
    end
end

params = struct('A',3.0,'B',0.9,'C',0.34, ...
    'G',1.1,'H',0.10,'nBoundary',160);
geom = make_boundary(params);
meshData = generate_mesh(geom,struct('targetH',0.11));
K = assemble_stiffness(meshData);
psi0 = initial_guess(meshData,K, ...
    struct('boundaryValue',0,'axisValue',-1));
profiles = make_profiles(struct( ...
    'pprime',struct('type','constant','value',0), ...
    'FFprime',struct('type','constant','value',-1)));
solverOptions = struct('tolerance',1e-11,'residualTolerance',1e-11, ...
    'maxIterations',10,'omega',1,'boundaryValue',0, ...
    'axisMode','min','verbose',false);
[psi,result] = solve_picard(meshData,profiles,psi0,solverOptions);
assert(result.converged,'Linear verification equilibrium did not converge.');

axisData = find_axis(meshData,psi, ...
    struct('mode','min','boundaryValue',0));
assert(inpolygon(axisData.R,axisData.Z,geom.R,geom.Z), ...
    'Refined magnetic axis lies outside the LCFS.');
assert(axisData.psi < 0,'Expected a negative interior flux axis.');

topology = check_topology(meshData,psi,axisData,struct( ...
    'normalizedLevels',(0.1:0.1:0.9).', ...
    'boundaryValue',0,'failOnFailure',false));
assert(topology.contourPassed, ...
    'Closed nested contour verification failed.');
assert(topology.nOPoints==1, ...
    'Expected exactly one discrete O-point candidate.');
assert(topology.nXPoints==0, ...
    'Expected no discrete X-point candidate.');
assert(topology.passed,'Complete topology verification failed.');

% For psi=R, the exact P1 gradient is dpsi/dR=1 and dpsi/dZ=0.
psiManufactured = meshData.nodes(:,1);
field = compute_B(meshData,psiManufactured,2.0);
exactElementBZ = 1./field.element.centroid(:,1);
exactNodeBZ = 1./meshData.nodes(:,1);
elementBRError = norm(field.element.BR,inf);
elementBZError = norm(field.element.BZ-exactElementBZ,inf);
nodeBZError = norm(field.node.BZ-exactNodeBZ,inf);
elementBphiError = norm(field.element.Bphi ...
    -2./field.element.centroid(:,1),inf);
assert(elementBRError<1e-12 && elementBZError<1e-12 ...
        && nodeBZError<1e-12 && elementBphiError<1e-12, ...
    'Manufactured magnetic-field verification failed.');

% Pure toroidal unit field: cylindrical geometry gives kappa=-e_R/R.
% Combine it with manufactured psi=R and p'=2 to verify
% kappa dot grad(p)=-2/R, including the cylindrical basis derivative.
toroidalField = field;
toroidalField.node.BR = zeros(meshData.nNodes,1);
toroidalField.node.BZ = zeros(meshData.nNodes,1);
toroidalField.node.Bp = zeros(meshData.nNodes,1);
toroidalField.node.Bphi = ones(meshData.nNodes,1);
toroidalField.node.Btotal = ones(meshData.nNodes,1);
toroidalField.element.BR = zeros(meshData.nElements,1);
toroidalField.element.BZ = zeros(meshData.nElements,1);
toroidalField.element.Bp = zeros(meshData.nElements,1);
toroidalField.element.Bphi = ones(meshData.nElements,1);
toroidalField.element.Btotal = ones(meshData.nElements,1);
pressureProfiles = make_profiles(struct( ...
    'pprime',struct('type','constant','value',2), ...
    'FFprime',struct('type','constant','value',0)));
curvature = compute_magnetic_curvature( ...
    meshData,psiManufactured,toroidalField,pressureProfiles,struct( ...
        'psiAxis',min(psiManufactured), ...
        'boundaryValue',max(psiManufactured)));
centroidR = field.element.centroid(:,1);
curvatureRError = norm(curvature.element.kappaR+1./centroidR,inf);
curvaturePhiError = norm(curvature.element.kappaPhi,inf);
curvatureZError = norm(curvature.element.kappaZ,inf);
goodCurvatureError = norm( ...
    curvature.element.kappaDotGradP+2./centroidR,inf);
assert(curvatureRError<1e-12 && curvaturePhiError<1e-12 && ...
        curvatureZError<1e-12 && goodCurvatureError<1e-12, ...
    'Manufactured magnetic-curvature verification failed.');

report = struct();
report.nNodes = meshData.nNodes;
report.nElements = meshData.nElements;
report.axisR = axisData.R;
report.axisZ = axisData.Z;
report.axisPsi = axisData.psi;
report.axisRefined = axisData.refinementAccepted;
report.nTopologyLevels = numel(topology.levelChecks);
report.nOPoints = topology.nOPoints;
report.nXPoints = topology.nXPoints;
report.topologyPassed = topology.passed;
report.elementBRError = elementBRError;
report.elementBZError = elementBZError;
report.nodeBZError = nodeBZError;
report.elementBphiError = elementBphiError;
report.curvatureRError = curvatureRError;
report.curvaturePhiError = curvaturePhiError;
report.curvatureZError = curvatureZError;
report.kappaDotGradPError = goodCurvatureError;

fprintf(['Post verification: axis=(%.6f, %.6f), psi_axis=%.6g, ' ...
         'levels=%d, O=%d, X=%d, BZ error=%.3e, curvature error=%.3e\n'], ...
    report.axisR,report.axisZ,report.axisPsi, ...
    report.nTopologyLevels,report.nOPoints,report.nXPoints, ...
    report.elementBZError,report.curvatureRError);
fprintf('post verification passed.\n');

if doPlot
    plot_flux(meshData,psi,struct('axisData',axisData, ...
        'boundaryValue',0,'title','Verified concave C-shape equilibrium'));
end

end
