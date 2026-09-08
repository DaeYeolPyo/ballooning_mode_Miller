function report = verify_physics(doPlot)
%VERIFY_PHYSICS Verify flux normalization, profiles, and FEM source coupling.

if nargin < 1
    doPlot = true;
end
validateattributes(doPlot, {'logical', 'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'doPlot');
doPlot = logical(doPlot);

flux = struct( ...
    'psiAxis', -2.0, ...
    'psiBoundary', 0.0, ...
    'rangePolicy', 'error', ...
    'tolerance', 1e-10);
psi = [-2.0, -1.0, 0.0];
[psiNormalized, normInfo] = normalize_flux( ...
    psi, flux.psiAxis, flux.psiBoundary, struct( ...
    'rangePolicy', flux.rangePolicy, 'tolerance', flux.tolerance));
normalizationError = norm(psiNormalized-[0, 0.5, 1], inf);
assert(normalizationError < 1e-15, ...
    'Normalized flux did not map axis/midpoint/boundary to 0/0.5/1.');

clipped = normalize_flux([-2.2, 0.2], -2, 0, ...
    struct('rangePolicy', 'clip'));
assert(isequal(clipped, [0, 1]), 'Flux clipping failed.');

constantP = pprime(psi, -4.0);
assert(isequal(constantP, -4*ones(size(psi))), ...
    'Numeric constant pprime failed.');

pShape = struct('type', 'normalized-power', ...
    'axisValue', -4.0, 'exponent', 2.0);
pShapeValue = pprime(psi, pShape, flux);
assert(norm(pShapeValue-[-4, -1, 0], inf) < 1e-14, ...
    'Normalized-power pprime failed.');

pPhysical = struct('type', 'pressure-power', ...
    'axisPressure', 100.0, 'boundaryPressure', 10.0, 'exponent', 2.0);
pPhysicalValue = pprime(psi, pPhysical, flux);
assert(norm(pPhysicalValue-[-90, -45, 0], inf) < 1e-14, ...
    'Analytic pressure-power derivative failed.');

FFShape = struct('type', 'normalized-power', ...
    'axisValue', 6.0, 'exponent', 1.0);
FFShapeValue = FFprime(psi, FFShape, flux);
assert(norm(FFShapeValue-[6, 3, 0], inf) < 1e-14, ...
    'Normalized-power FFprime failed.');

FPhysical = struct('type', 'F-power', ...
    'axisF', 3.0, 'boundaryF', 2.0, 'exponent', 2.0);
FPhysicalValue = FFprime(psi, FPhysical, flux);
assert(norm(FPhysicalValue-[-3, -1.125, 0], inf) < 1e-14, ...
    'Analytic F-power FFprime failed.');

tabulatedP = struct('type','tabulated', ...
    'psiN',[0;0.5;1],'values',[0;-4;0],'method','pchip');
tabulatedPValue = pprime(psi,tabulatedP,flux);
assert(norm(tabulatedPValue-[0,-4,0],inf)<1e-14, ...
    'Tabulated pprime interpolation failed.');

tabulatedFF = struct('type','tabulated', ...
    'psiN',[0;0.5;1],'values',[-2;1;0],'method','pchip', ...
    'FValues',[3;2.5;2]);
tabulatedFFValue = FFprime(psi,tabulatedFF,flux);
assert(norm(tabulatedFFValue-[-2,1,0],inf)<1e-14, ...
    'Tabulated FFprime interpolation failed.');
tabulatedFValue = reconstruct_F(psi,2,tabulatedFF,flux);
assert(norm(tabulatedFValue-[3,2.5,2],inf)<1e-14, ...
    'Tabulated F reconstruction failed.');

overrides = struct();
overrides.flux = flux;
overrides.pprime = struct('type', 'constant', 'value', -2.5e5);
overrides.FFprime = struct('type', 'constant', 'value', 0.4);
profiles = make_profiles(overrides);

R = [2.0, 3.0, 4.0];
Z = zeros(size(R));
q = gs_source(R, Z, psi, profiles);
qExpected = profiles.mu0*R*overrides.pprime.value ...
          + overrides.FFprime.value./R;
sourceFormulaError = norm(q-qExpected, inf);
assert(sourceFormulaError < 1e-14, ...
    'gs_source does not match the divided weak-source formula.');

femCouplingError = NaN;
if exist('make_boundary', 'file') == 2 && ...
        exist('generate_mesh', 'file') == 2 && ...
        exist('assemble_rhs', 'file') == 2
    geom = make_boundary(struct( ...
        'A', 3.0, 'B', 0.9, 'C', 0.34, ...
        'G', 1.1, 'H', 0.10, 'nBoundary', 120));
    meshData = generate_mesh(geom, struct('targetH', 0.14));
    psiNodal = -1.5*ones(meshData.nNodes, 1);
    rhsPhysics = assemble_rhs(meshData, ...
        @(Rq,Zq,psiQ) gs_source(Rq,Zq,psiQ,profiles), psiNodal);
    rhsExplicit = assemble_rhs(meshData, ...
        @(Rq,~,~) profiles.mu0*Rq*overrides.pprime.value ...
                  + overrides.FFprime.value./Rq, psiNodal);
    femCouplingError = norm(rhsPhysics-rhsExplicit, inf) ...
                     / max(norm(rhsExplicit, inf), eps);
    assert(femCouplingError < 1e-14, ...
        'Physics callback and explicit FEM source assembly disagree.');
end

report = struct();
report.normalizationError = normalizationError;
report.normalizedRange = normInfo.rawRange;
report.sourceFormulaError = sourceFormulaError;
report.femCouplingError = femCouplingError;
report.pprimeAtAxis = pShapeValue(1);
report.FFprimeAtAxis = FFShapeValue(1);

fprintf(['Physics verification: normalization=%.3e, source=%.3e, ' ...
         'FEM coupling=%.3e\n'], ...
    report.normalizationError, report.sourceFormulaError, ...
    report.femCouplingError);
fprintf('physics verification passed.\n');

if doPlot
    psiGrid = linspace(flux.psiAxis, flux.psiBoundary, 201);
    psiNGrid = normalize_flux( ...
        psiGrid, flux.psiAxis, flux.psiBoundary, struct( ...
        'rangePolicy', flux.rangePolicy, 'tolerance', flux.tolerance));
    figure('Name', 'Grad-Shafranov source profiles', 'Color', 'w');
    subplot(1,2,1)
    plot(psiNGrid, pprime(psiGrid, pShape, flux), ...
        'LineWidth', 1.8)
    grid on
    xlabel('\psi_N')
    ylabel('p''(\psi)')
    title('Normalized-power p'' profile')

    subplot(1,2,2)
    plot(psiNGrid, FFprime(psiGrid, FFShape, flux), ...
        'LineWidth', 1.8)
    grid on
    xlabel('\psi_N')
    ylabel('FF''(\psi)')
    title('Normalized-power FF'' profile')
end

end
