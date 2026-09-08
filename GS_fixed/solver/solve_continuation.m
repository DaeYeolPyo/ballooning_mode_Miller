function continuation = solve_continuation( ...
    startParams, targetParams, lambdas, profiles, opts)
%SOLVE_CONTINUATION Track a fixed-boundary solution into C-shaping.
%
%   CONT = SOLVE_CONTINUATION(STARTPARAMS, TARGETPARAMS, LAMBDAS, PROFILES)
%   linearly interpolates the active make_boundary parameters from lambda=0
%   to lambda=1. Harmonic boundaries interpolate A,B,C,G,H; localized-c
%   boundaries interpolate A,B,D,beta,G,H; pointed boomerang-c boundaries
%   interpolate A,Rback,Rtip,Rend,Zend,pInner,pOuter. Each converged
%   Third-harmonic cosine-c boundaries interpolate R0,A,B,C,kappa. Each
%   Miller-D boundaries interpolate R0,a,kappa,delta. Each
%   Tabulated boundaries interpolate corresponding R and Z point arrays.
%   converged solution is transferred to the next mesh and used as its
%   Picard initial iterate.
%
%   LAMBDAS must be a strictly increasing vector in [0,1]. STARTPARAMS and
%   TARGETPARAMS are complete make_boundary parameter structs.
%
%   OPTS supports:
%       mesh         options passed to generate_mesh
%       initialGuess options passed to initial_guess
%       picard       options passed to solve_picard
%       verbose      print continuation diagnostics (default false)
%
%   CONT.steps(k) stores lambda, geom, mesh, psi, Picard result, and transfer
%   diagnostics. Keeping complete step data is intentional for branch and
%   topology inspection in the upcoming post-processing stage.

narginchk(4, 5);

require_function('make_boundary');
require_function('generate_mesh');
require_function('assemble_stiffness');
require_function('initial_guess');
require_function('solve_picard');
require_function('transfer_solution');

if nargin < 5 || isempty(opts)
    opts = struct();
end
opts = validate_options(opts);

validateattributes(lambdas, {'numeric'}, ...
    {'real', 'finite', 'vector', 'nonempty', '>=', 0, '<=', 1}, ...
    mfilename, 'lambdas');
lambdas = lambdas(:);
if any(diff(lambdas) <= 0)
    error('GS:solver:InvalidContinuationSequence', ...
        'lambdas must be strictly increasing.');
end

[startParams, startType] = canonical_boundary_params( ...
    startParams, 'startParams');
[targetParams, targetType] = canonical_boundary_params( ...
    targetParams, 'targetParams');
if ~strcmp(startType, targetType)
    error('GS:solver:BoundaryTypeMismatch', ...
        'startParams and targetParams must use the same boundary type.');
end
if strcmp(startType,'tabulated') && ...
        numel(startParams.R)~=numel(targetParams.R)
    error('GS:solver:TabulatedBoundarySizeMismatch', ...
        ['Tabulated start and target boundaries must contain the same ' ...
         'number of corresponding points.']);
end

nSteps = numel(lambdas);
steps = repmat(struct( ...
    'lambda', [], 'params', [], 'geom', [], 'mesh', [], ...
    'psi', [], 'solver', [], 'transfer', []), nSteps, 1);

previousMesh = [];
previousPsi = [];
for k = 1:nSteps
    lambda = lambdas(k);
    params = interpolate_params( ...
        startParams, targetParams, startType, lambda);
    geom = make_boundary(params);
    meshData = generate_mesh(geom, opts.mesh);
    K = assemble_stiffness(meshData);

    if k == 1
        [psi0, initialInfo] = initial_guess( ...
            meshData, K, opts.initialGuess);
        transferInfo = struct( ...
            'method', 'elliptic initial guess', ...
            'nInterpolated', 0, 'nExtrapolated', 0, ...
            'extrapolatedFraction', 0, ...
            'initialGuess', initialInfo);
    else
        boundaryValue = picard_boundary_value(opts.picard);
        [psi0, transferInfo] = transfer_solution( ...
            previousMesh, previousPsi, meshData, boundaryValue);
        transferInfo.method = 'P1 barycentric plus nearest extrapolation';
    end

    [psi, solverResult] = solve_picard( ...
        meshData, profiles, psi0, opts.picard);

    steps(k).lambda = lambda;
    steps(k).params = params;
    steps(k).geom = geom;
    steps(k).mesh = meshData;
    steps(k).psi = psi;
    steps(k).solver = solverResult;
    steps(k).transfer = transferInfo;

    if opts.verbose
        fprintf(['Continuation %d/%d: lambda=%.4f, nodes=%d, ' ...
                 'Picard=%d, residual=%.3e, psi_axis=%.6g\n'], ...
            k, nSteps, lambda, meshData.nNodes, ...
            solverResult.nIterations, solverResult.finalResidual, ...
            solverResult.axisValue);
    end

    previousMesh = meshData;
    previousPsi = psi;
end

continuation = struct();
continuation.lambdas = lambdas;
continuation.steps = steps;
continuation.finalGeometry = steps(end).geom;
continuation.finalMesh = steps(end).mesh;
continuation.finalPsi = steps(end).psi;
continuation.finalSolver = steps(end).solver;
continuation.options = opts;

end


function params = interpolate_params( ...
        startParams, targetParams, type, lambda)
if strcmp(type, 'harmonic')
    shapeFields = {'A', 'B', 'C', 'G', 'H'};
elseif strcmp(type, 'localized-c')
    shapeFields = {'A', 'B', 'D', 'beta', 'G', 'H'};
elseif strcmp(type, 'cosine-c')
    shapeFields = {'R0', 'A', 'B', 'C', 'kappa'};
elseif strcmp(type,'miller-d')
    shapeFields = {'R0','a','kappa','delta'};
elseif strcmp(type,'tabulated')
    shapeFields = {'R','Z'};
else
    shapeFields = {'A', 'Rback', 'Rtip', 'Rend', 'Zend', ...
        'pInner', 'pOuter'};
end
params = startParams;
for k = 1:numel(shapeFields)
    name = shapeFields{k};
    params.(name) = (1-lambda)*startParams.(name) ...
                  + lambda*targetParams.(name);
end
params.nBoundary = round((1-lambda)*startParams.nBoundary ...
                       + lambda*targetParams.nBoundary);
if strcmp(type, 'boomerang-c')
    params.nBoundary = 2*round(params.nBoundary/2);
end
end


function [params, type] = canonical_boundary_params(params, name)
if ~isstruct(params) || ~isscalar(params)
    error('GS:solver:InvalidBoundaryParameters', ...
        '%s must be a scalar boundary-parameter struct.', name);
end
if isfield(params, 'type')
    type = normalize_boundary_type(params.type);
else
    type = 'harmonic';
end
params.type = type;
if strcmp(type, 'harmonic')
    required = {'type', 'A', 'B', 'C', 'G', 'H', 'nBoundary'};
elseif strcmp(type, 'localized-c')
    required = {'type', 'A', 'B', 'D', 'beta', 'G', 'H', 'nBoundary'};
elseif strcmp(type, 'cosine-c')
    required = {'type', 'R0', 'A', 'B', 'C', 'kappa', 'nBoundary'};
elseif strcmp(type,'miller-d')
    required = {'type','R0','a','kappa','delta','nBoundary'};
elseif strcmp(type,'tabulated')
    required = {'type','R','Z','nBoundary'};
else
    required = {'type', 'A', 'Rback', 'Rtip', 'Rend', 'Zend', ...
        'pInner', 'pOuter', 'nBoundary'};
end
if ~all(isfield(params, required)) || ...
        ~isempty(setdiff(fieldnames(params), required))
    error('GS:solver:InvalidBoundaryParameters', ...
        '%s does not contain exactly the fields required for type %s.', ...
        name, type);
end
if strcmp(type,'tabulated') && numel(params.R)~=numel(params.Z)
    error('GS:solver:InvalidBoundaryParameters', ...
        '%s tabulated R and Z arrays must have equal lengths.',name);
end
end


function type = normalize_boundary_type(value)
if isstring(value) && isscalar(value)
    value = char(value);
end
if ~ischar(value) || ~isrow(value)
    error('GS:solver:InvalidBoundaryType', ...
        ['Boundary type must be harmonic, localized-c, boomerang-c, ' ...
         'cosine-c, miller-d, or tabulated.']);
end
type = lower(strrep(strtrim(value), '_', '-'));
if ismember(type, {'localized', 'local-c'})
    type = 'localized-c';
end
if ismember(type, {'boomerang', 'crescent'})
    type = 'boomerang-c';
end
if ismember(type, {'cosine', 'fourier-c', 'harmonic3', 'third-harmonic'})
    type = 'cosine-c';
end
if ismember(type,{'miller','d-shape','dshape','miller-dshape'})
    type = 'miller-d';
end
if ismember(type,{'points','geqdsk','sampled'})
    type = 'tabulated';
end
if ~ismember(type, ...
        {'harmonic','localized-c','boomerang-c','cosine-c','miller-d', ...
         'tabulated'})
    error('GS:solver:InvalidBoundaryType', ...
        ['Boundary type must be harmonic, localized-c, boomerang-c, ' ...
         'cosine-c, miller-d, or tabulated.']);
end
end


function value = picard_boundary_value(picardOpts)
if isfield(picardOpts, 'boundaryValue')
    value = picardOpts.boundaryValue;
else
    value = 0.0;
end
end


function opts = validate_options(opts)
if ~isstruct(opts) || ~isscalar(opts)
    error('GS:solver:InvalidContinuationOptions', ...
        'opts must be a scalar struct.');
end
defaults = struct( ...
    'mesh', struct(), ...
    'initialGuess', struct(), ...
    'picard', struct(), ...
    'verbose', false);
unknown = setdiff(fieldnames(opts), fieldnames(defaults));
if ~isempty(unknown)
    error('GS:solver:UnknownContinuationOption', ...
        'Unknown continuation option: %s', strjoin(unknown, ', '));
end
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(opts, names{k})
        opts.(names{k}) = defaults.(names{k});
    end
end
for name = {'mesh', 'initialGuess', 'picard'}
    value = opts.(name{1});
    if ~isstruct(value) || ~isscalar(value)
        error('GS:solver:InvalidContinuationSuboptions', ...
            'opts.%s must be a scalar struct.', name{1});
    end
end
validateattributes(opts.verbose, {'logical', 'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, 'opts.verbose');
opts.verbose = logical(opts.verbose);

boundaryValue = picard_boundary_value(opts.picard);
if isfield(opts.initialGuess, 'boundaryValue') && ...
        opts.initialGuess.boundaryValue ~= boundaryValue
    error('GS:solver:InconsistentBoundaryValue', ...
        ['initialGuess.boundaryValue and picard.boundaryValue must agree ' ...
         'during continuation.']);
end
opts.initialGuess.boundaryValue = boundaryValue;
end


function require_function(name)
if exist(name, 'file') ~= 2
    error('GS:solver:MissingDependency', ...
        'Required function %s is not on the MATLAB path.', name);
end
end
