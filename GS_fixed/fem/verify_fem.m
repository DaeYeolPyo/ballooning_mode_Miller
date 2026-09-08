function report = verify_fem(doPlot)
%VERIFY_FEM Verify P1 assembly and fixed-boundary elimination.
%
%   REPORT = VERIFY_FEM() runs algebraic checks on a concave C-shaped mesh,
%   solves the manufactured problem psi=R with q=1/R^2, and plots numerical
%   and nodal-error fields. VERIFY_FEM(false) suppresses the plot.

if nargin < 1
    doPlot = true;
end
validateattributes(doPlot, {'logical', 'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'doPlot');
doPlot = logical(doPlot);

if exist('make_boundary', 'file') ~= 2 || ...
        exist('generate_mesh', 'file') ~= 2
    error('GS:fem:GeometryNotOnPath', ...
        ['Add both geometry/ and fem/ to the MATLAB path before running ' ...
         'verify_fem.']);
end

params = struct( ...
    'A', 3.0, 'B', 0.9, 'C', 0.34, ...
    'G', 1.1, 'H', 0.10, 'nBoundary', 160);
geom = make_boundary(params);
meshData = generate_mesh(geom, struct('targetH', 0.10));

[K, stiffnessInfo] = assemble_stiffness(meshData);
nNodes = meshData.nNodes;

symmetryError = norm(K-K.', 'fro')/norm(K, 'fro');
nullResidual = norm(K*ones(nNodes,1), inf)/norm(K, 'fro');
assert(symmetryError < 1e-14, 'The stiffness matrix is not symmetric.');
assert(nullResidual < 1e-13, ...
    'The constant-field null residual is too large.');

[constantLoad, loadInfo] = assemble_rhs(meshData, 1.0);
nodalLoad = assemble_rhs(meshData, ones(nNodes,1));
callbackLoad = assemble_rhs(meshData, @(R,~) ones(size(R)));
psiProbe = meshData.nodes(:,1) + 0.2*meshData.nodes(:,2);
nodalProbeLoad = assemble_rhs(meshData, psiProbe);
nonlinearProbeLoad = assemble_rhs( ...
    meshData, @(~,~,psi) psi, psiProbe);
loadAreaError = abs(sum(constantLoad)-geom.area)/geom.area;
assert(loadAreaError < 1e-13, ...
    'A unit source must integrate to the LCFS polygon area.');
assert(norm(constantLoad-nodalLoad, inf) < 1e-14, ...
    'Scalar and nodal unit sources must assemble identically.');
assert(norm(constantLoad-callbackLoad, inf) < 1e-14, ...
    'Scalar and callback unit sources must assemble identically.');
nonlinearCallbackError = norm(nodalProbeLoad-nonlinearProbeLoad, inf);
assert(nonlinearCallbackError < 1e-14, ...
    'Nodal and nonlinear-callback sources must assemble identically.');

% Constant Dirichlet data with zero source must reproduce the same constant
% throughout the domain exactly up to the linear solve tolerance.
[Kconstant, rhsConstant] = apply_dirichlet( ...
    meshData, K, zeros(nNodes,1), 2.5);
[~, cholStatus] = chol(Kconstant);
assert(cholStatus == 0, ...
    'The Dirichlet-modified stiffness matrix must be positive definite.');
psiConstant = Kconstant\rhsConstant;
constantSolutionError = max(abs(psiConstant-2.5));
assert(constantSolutionError < 1e-11, ...
    'Constant Dirichlet data were not reproduced accurately.');

% Manufactured strong/weak problem:
%   psi_exact = R,
%   -div((1/R) grad(psi_exact)) = 1/R^2.
rhsManufactured = assemble_rhs(meshData, @(R,~) 1./R.^2);
[Kmanufactured, rhsManufacturedBc, bc] = apply_dirichlet( ...
    meshData, K, rhsManufactured, @(R,~) R);
psiNumerical = Kmanufactured\rhsManufacturedBc;
psiExact = meshData.nodes(:,1);
errorField = psiNumerical-psiExact;
relativeL2Error = norm(errorField)/norm(psiExact);
relativeInfinityError = norm(errorField, inf)/norm(psiExact, inf);
freeResidual = norm(K(bc.freeNodes,:)*psiNumerical ...
    - rhsManufactured(bc.freeNodes), inf);

assert(relativeL2Error < 2e-4, ...
    'The manufactured psi=R solution error is unexpectedly large.');
assert(freeResidual < 1e-10, ...
    'The manufactured solution has an unexpectedly large free-node residual.');

report = struct();
report.nNodes = meshData.nNodes;
report.nElements = meshData.nElements;
report.nNonzeros = stiffnessInfo.nNonzeros;
report.symmetryError = symmetryError;
report.constantNullResidual = nullResidual;
report.loadAreaError = loadAreaError;
report.loadIntegral = loadInfo.integratedSource;
report.nonlinearCallbackError = nonlinearCallbackError;
report.constantSolutionError = constantSolutionError;
report.manufacturedRelativeL2Error = relativeL2Error;
report.manufacturedRelativeInfinityError = relativeInfinityError;
report.manufacturedFreeResidual = freeResidual;

fprintf(['FEM verification: nodes=%d, elements=%d, nnz(K)=%d, ' ...
         'sym=%.3e, null=%.3e, area=%.3e, manufactured L2=%.3e\n'], ...
    report.nNodes, report.nElements, report.nNonzeros, ...
    report.symmetryError, report.constantNullResidual, ...
    report.loadAreaError, report.manufacturedRelativeL2Error);
fprintf('fem verification passed.\n');

if doPlot
    figure('Name', 'P1 FEM manufactured-solution verification', ...
        'Color', 'w');
    subplot(1,2,1)
    trisurf(meshData.elements, meshData.nodes(:,1), ...
        meshData.nodes(:,2), psiNumerical, 'EdgeColor', 'none');
    view(2)
    axis equal tight
    colorbar
    xlabel('R')
    ylabel('Z')
    title('Numerical \psi for exact \psi=R')

    subplot(1,2,2)
    trisurf(meshData.elements, meshData.nodes(:,1), ...
        meshData.nodes(:,2), errorField, 'EdgeColor', 'none');
    view(2)
    axis equal tight
    colorbar
    xlabel('R')
    ylabel('Z')
    title('Nodal error')
end

end
