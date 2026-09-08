function [psi, result] = solve_picard(meshData, profiles, psi0, opts)
%SOLVE_PICARD Solve the nonlinear fixed-boundary GS equation by Picard steps.
%
%   [PSI, RESULT] = SOLVE_PICARD(MESHDATA, PROFILES, PSI0, OPTS) repeatedly
%   evaluates
%
%       q(psi^n) = mu0*R*p'(psi^n) + FF'(psi^n)/R,
%
%   solves the P1 FEM system, and applies under-relaxation. The LCFS flux is
%   a single scalar because a fixed plasma boundary is a magnetic surface.
%
%   Supported OPTS fields:
%       tolerance          relative iterate tolerance       (1e-8)
%       residualTolerance  relative free-node residual      (1e-8)
%       maxIterations      maximum Picard steps             (300)
%       omega              under-relaxation in (0,1]        (0.5)
%       boundaryValue      constant LCFS psi                (0)
%       axisMode           'auto', 'min', or 'max'          ('auto')
%       updateFluxAxis     update profile psiAxis each step (true)
%       failOnNonconvergence                                (true)
%       verbose            print iteration diagnostics      (false)
%       printEvery         verbose output interval          (10)
%
%   RESULT contains convergence histories, final source residual, estimated
%   nodal magnetic-axis value/location, and the stiffness matrix.

narginchk(3, 4);

require_function('assemble_stiffness');
require_function('assemble_rhs');
require_function('apply_dirichlet');
require_function('gs_source');

if nargin < 4 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts);

requiredMesh = {'nodes', 'elements', 'boundaryNodes', 'interiorNodes'};
if ~isstruct(meshData) || ~isscalar(meshData) || ...
        ~all(isfield(meshData, requiredMesh))
    error('GS:solver:InvalidMesh', ...
        'meshData must be a scalar mesh returned by generate_mesh.');
end
nNodes = size(meshData.nodes, 1);
if isempty(meshData.interiorNodes)
    error('GS:solver:NoInteriorNodes', ...
        'The fixed-boundary mesh must contain interior nodes.');
end
validateattributes(psi0, {'numeric'}, ...
    {'real', 'finite', 'vector', 'numel', nNodes}, mfilename, 'psi0');
psi = psi0(:);
psi(meshData.boundaryNodes) = opts.boundaryValue;

if ~isstruct(profiles) || ~isscalar(profiles) || ...
        ~all(isfield(profiles, {'mu0', 'pprime', 'FFprime'}))
    error('GS:solver:InvalidProfiles', ...
        'profiles must be a scalar configuration from make_profiles.');
end

[K, stiffnessInfo] = assemble_stiffness(meshData);
[Kbc, ~, bc] = apply_dirichlet( ...
    meshData, K, zeros(nNodes,1), opts.boundaryValue);
linearSolver = decomposition(Kbc, 'chol');

updateHistory = NaN(opts.maxIterations, 1);
residualHistory = NaN(opts.maxIterations, 1);
axisHistory = NaN(opts.maxIterations, 1);
sourceNormHistory = NaN(opts.maxIterations, 1);
converged = false;

for iteration = 1:opts.maxIterations
    [profilesIteration, axisValue] = iteration_profiles( ...
        profiles, psi, meshData.interiorNodes, opts);
    rhs = assemble_rhs(meshData, ...
        @(R,Z,psiQ) gs_source(R,Z,psiQ,profilesIteration), psi);
    [~, rhsBc] = apply_dirichlet( ...
        meshData, K, rhs, opts.boundaryValue);
    psiLinear = linearSolver\rhsBc;

    psiNew = psi;
    psiNew(bc.freeNodes) = (1-opts.omega)*psi(bc.freeNodes) ...
                         + opts.omega*psiLinear(bc.freeNodes);
    psiNew(bc.boundaryNodes) = opts.boundaryValue;

    [profilesNew, axisValueNew] = iteration_profiles( ...
        profiles, psiNew, meshData.interiorNodes, opts);
    rhsNew = assemble_rhs(meshData, ...
        @(R,Z,psiQ) gs_source(R,Z,psiQ,profilesNew), psiNew);
    residual = K(bc.freeNodes,:)*psiNew-rhsNew(bc.freeNodes);

    updateScale = max(norm(psiNew(bc.freeNodes), 2), eps);
    updateHistory(iteration) = ...
        norm(psiNew(bc.freeNodes)-psi(bc.freeNodes), 2)/updateScale;
    residualScale = max([norm(K(bc.freeNodes,:)*psiNew, 2), ...
                         norm(rhsNew(bc.freeNodes), 2), eps]);
    residualHistory(iteration) = norm(residual, 2)/residualScale;
    axisHistory(iteration) = axisValueNew;
    sourceNormHistory(iteration) = norm(rhsNew(bc.freeNodes), 2);

    if opts.verbose && (iteration == 1 || ...
            mod(iteration, opts.printEvery) == 0)
        fprintf('Picard %4d: update=%.3e, residual=%.3e, psi_axis=%.6g\n', ...
            iteration, updateHistory(iteration), ...
            residualHistory(iteration), axisValueNew);
    end

    psi = psiNew;
    if updateHistory(iteration) <= opts.tolerance && ...
            residualHistory(iteration) <= opts.residualTolerance
        converged = true;
        break
    end

    if ~isfinite(axisValue) || ~isfinite(axisValueNew)
        error('GS:solver:NonFiniteAxisEstimate', ...
            'The Picard iteration produced a non-finite axis estimate.');
    end
end

nIterations = iteration;
updateHistory = updateHistory(1:nIterations);
residualHistory = residualHistory(1:nIterations);
axisHistory = axisHistory(1:nIterations);
sourceNormHistory = sourceNormHistory(1:nIterations);

[finalAxisValue, finalAxisNode, selectedAxisMode] = estimate_axis( ...
    psi, meshData.interiorNodes, opts.boundaryValue, opts.axisMode);

result = struct();
result.converged = converged;
result.nIterations = nIterations;
result.updateHistory = updateHistory;
result.residualHistory = residualHistory;
result.axisHistory = axisHistory;
result.sourceNormHistory = sourceNormHistory;
result.finalUpdate = updateHistory(end);
result.finalResidual = residualHistory(end);
result.axisMode = selectedAxisMode;
result.axisValue = finalAxisValue;
result.axisNode = finalAxisNode;
result.axisPosition = meshData.nodes(finalAxisNode,:);
result.boundaryValue = opts.boundaryValue;
result.stiffness = K;
result.stiffnessInfo = stiffnessInfo;
result.options = opts;

if opts.verbose
    fprintf(['Picard finished: converged=%d, iterations=%d, ' ...
             'update=%.3e, residual=%.3e\n'], ...
        converged, nIterations, result.finalUpdate, result.finalResidual);
end

if ~converged && opts.failOnNonconvergence
    error('GS:solver:PicardDidNotConverge', ...
        ['Picard iteration did not converge in %d steps. Final update is ' ...
         '%.3e and residual is %.3e. Reduce omega or use continuation.'], ...
        nIterations, result.finalUpdate, result.finalResidual);
end

end


function [profilesIteration, axisValue] = iteration_profiles( ...
    profiles, psi, interiorNodes, opts)
[axisValue, ~] = estimate_axis( ...
    psi, interiorNodes, opts.boundaryValue, opts.axisMode);
axisScale = max([1, abs(axisValue), abs(opts.boundaryValue)]);
if abs(axisValue-opts.boundaryValue) <= 100*eps(axisScale)
    error('GS:solver:CollapsedFluxRange', ...
        ['The interior flux extremum has collapsed onto the boundary value. ' ...
         'Use a nonconstant initial guess with the expected source sign.']);
end

profilesIteration = profiles;
if opts.updateFluxAxis
    if ~isfield(profilesIteration, 'flux') || ...
            ~isstruct(profilesIteration.flux)
        profilesIteration.flux = struct();
    end
    profilesIteration.flux.psiAxis = axisValue;
    profilesIteration.flux.psiBoundary = opts.boundaryValue;
    if ~isfield(profilesIteration.flux, 'rangePolicy')
        profilesIteration.flux.rangePolicy = 'error';
    end
    if ~isfield(profilesIteration.flux, 'tolerance')
        profilesIteration.flux.tolerance = 1e-10;
    end
end
end


function [axisValue, axisNode, selectedMode] = estimate_axis( ...
    psi, interiorNodes, boundaryValue, mode)
interiorValues = psi(interiorNodes);
[minimum, minimumIndex] = min(interiorValues);
[maximum, maximumIndex] = max(interiorValues);

if strcmp(mode, 'auto')
    if abs(minimum-boundaryValue) >= abs(maximum-boundaryValue)
        selectedMode = 'min';
    else
        selectedMode = 'max';
    end
else
    selectedMode = mode;
end

if strcmp(selectedMode, 'min')
    axisValue = minimum;
    axisNode = interiorNodes(minimumIndex);
else
    axisValue = maximum;
    axisNode = interiorNodes(maximumIndex);
end
end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:solver:InvalidPicardOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct( ...
    'tolerance', 1e-8, ...
    'residualTolerance', 1e-8, ...
    'maxIterations', 300, ...
    'omega', 0.5, ...
    'boundaryValue', 0.0, ...
    'axisMode', 'auto', ...
    'updateFluxAxis', true, ...
    'failOnNonconvergence', true, ...
    'verbose', false, ...
    'printEvery', 10);
unknown = setdiff(fieldnames(opts), fieldnames(defaults));
if ~isempty(unknown)
    error('GS:solver:UnknownPicardOption', ...
        'Unknown Picard option: %s', strjoin(unknown, ', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts, names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end

validateattributes(opts.tolerance, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, 'opts.tolerance');
validateattributes(opts.residualTolerance, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, ...
    mfilename, 'opts.residualTolerance');
validateattributes(opts.maxIterations, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
    mfilename, 'opts.maxIterations');
validateattributes(opts.omega, {'numeric'}, ...
    {'real', 'finite', 'scalar', '>', 0, '<=', 1}, ...
    mfilename, 'opts.omega');
validateattributes(opts.boundaryValue, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'opts.boundaryValue');
validateattributes(opts.printEvery, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
    mfilename, 'opts.printEvery');
opts.updateFluxAxis = logical_scalar(opts.updateFluxAxis, 'updateFluxAxis');
opts.failOnNonconvergence = logical_scalar( ...
    opts.failOnNonconvergence, 'failOnNonconvergence');
opts.verbose = logical_scalar(opts.verbose, 'verbose');
opts.axisMode = axis_mode(opts.axisMode);
end


function value = logical_scalar(value, name)
validateattributes(value, {'logical', 'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, ['opts.' name]);
value = logical(value);
end


function mode = axis_mode(value)
if isstring(value)
    if ~isscalar(value)
        error('GS:solver:InvalidAxisMode', ...
            'axisMode must be auto, min, or max.');
    end
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:solver:InvalidAxisMode', ...
        'axisMode must be auto, min, or max.');
end
mode = lower(strtrim(value));
if ~ismember(mode, {'auto', 'min', 'max'})
    error('GS:solver:InvalidAxisMode', ...
        'axisMode must be auto, min, or max.');
end
end


function require_function(name)
if exist(name, 'file') ~= 2
    error('GS:solver:MissingDependency', ...
        'Required function %s is not on the MATLAB path.', name);
end
end
