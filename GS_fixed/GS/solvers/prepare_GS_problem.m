function problem = prepare_GS_problem(input)
%[PREPARE_GS_PROBLEM]
%
% Construct all profile-independent objects required by the
% fixed-boundary Grad-Shafranov solver.
%
% This function prepares:
%
%   1. P1 and P3 meshes
%   2. P3 quadrature
%   3. Global stiffness matrix
%   4. Dirichlet node partition
%   5. Reduced-matrix Cholesky factorization
%   6. Default initial psi
%   7. Picard options
%   8. Common profile and q-evaluation grids
%
% The returned objects can be reused during every outer iteration.
%
%
% INPUT
%
%   input
%       Existing solve_GS input structure.
%
%
% OUTPUT
%
%   problem
%       Prepared problem structure containing:
%
%       problem.mesh
%       problem.meshP1
%       problem.quad
%       problem.operator
%       problem.boundary
%       problem.dim
%       problem.grid
%       problem.profile
%       problem.initial
%       problem.options
%       problem.qEvaluation
%       problem.output
%       problem.timing

    narginchk(1, 1);

    if ~isstruct(input) || ~isscalar(input)
        error('GS:PrepareProblem:InputType', ...
            'input must be a scalar structure.');
    end

    %==============================================================
    % Required top-level structures
    %==============================================================
    assert_required_fields(input, ...
        {'grid', 'dim', 'boundary', 'Picard', 'output'}, ...
        'input');

    assert_required_fields(input.grid, ...
        {'nBoundary', 'NR', 'NZ', ...
         'quadratureDegree', 'nProfile'}, ...
        'input.grid');

    assert_required_fields(input.dim, ...
        {'R0', 'B0', 'aMinor'}, ...
        'input.dim');

    assert_required_fields(input.boundary, ...
        {'R', 'Z', 'psiBoundary'}, ...
        'input.boundary');

    assert_required_fields(input.Picard, ...
        {'omega', ...
         'maxIterations', ...
         'updateTolerance', ...
         'residualTolerance', ...
         'axisMode', ...
         'symmetryTolerance', ...
         'fluxSpanTolerance'}, ...
        'input.Picard');

    assert_required_fields(input.output, ...
        {'verbose', 'checkTime', 'showPlot'}, ...
        'input.output');

    %==============================================================
    % Grid inputs
    %==============================================================
    nBoundary = input.grid.nBoundary;
    NR = input.grid.NR;
    NZ = input.grid.NZ;

    quadratureDegree = ...
        input.grid.quadratureDegree;

    nProfile = input.grid.nProfile;

    validateattributes(nBoundary, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 3});

    validateattributes(NR, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2});

    validateattributes(NZ, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2});

    validateattributes(quadratureDegree, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', ...
         '>=', 1, '<=', 20});

    validateattributes(nProfile, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2});

    %==============================================================
    % Dimensional inputs
    %==============================================================
    R0 = input.dim.R0;
    B0 = input.dim.B0;
    aMinor = input.dim.aMinor;

    validateattributes(R0, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(B0, {'numeric'}, ...
        {'real', 'finite', 'scalar'});

    validateattributes(aMinor, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    %==============================================================
    % Boundary inputs
    %==============================================================
    R_bnd = input.boundary.R(:);
    Z_bnd = input.boundary.Z(:);

    psiBoundary = ...
        input.boundary.psiBoundary;

    if numel(R_bnd) ~= numel(Z_bnd)
        error('GS:PrepareProblem:BoundarySize', ...
            'Boundary R and Z arrays must have equal lengths.');
    end

    if numel(R_bnd) ~= nBoundary
        error('GS:PrepareProblem:BoundaryCount', ...
            ['input.grid.nBoundary does not match the number ', ...
             'of supplied boundary points.']);
    end

    if any(~isfinite(R_bnd)) || any(~isfinite(Z_bnd))
        error('GS:PrepareProblem:InvalidBoundary', ...
            'Boundary coordinates contain NaN or Inf.');
    end

    if any(R_bnd <= 0)
        error('GS:PrepareProblem:BoundaryRadius', ...
            'All boundary radii must be positive.');
    end

    validateattributes(psiBoundary, {'numeric'}, ...
        {'real', 'finite', 'scalar'});

    boundaryArea = polyarea(R_bnd, Z_bnd);

    coordinateScale = max([ ...
        max(R_bnd)-min(R_bnd), ...
        max(Z_bnd)-min(Z_bnd), ...
        1]);

    if ~isfinite(boundaryArea) || ...
            boundaryArea <= 100*eps(coordinateScale^2)
        error('GS:PrepareProblem:BoundaryArea', ...
            'The prescribed plasma boundary has negligible area.');
    end

    %==============================================================
    % Output options
    %==============================================================
    validateattributes(input.output.verbose, ...
        {'logical', 'numeric'}, ...
        {'real', 'finite', 'scalar'});

    validateattributes(input.output.checkTime, ...
        {'logical', 'numeric'}, ...
        {'real', 'finite', 'scalar'});

    validateattributes(input.output.showPlot, ...
        {'logical', 'numeric'}, ...
        {'real', 'finite', 'scalar'});

    verbose = logical(input.output.verbose);
    checkTime = logical(input.output.checkTime);
    showPlot = logical(input.output.showPlot);

    %==============================================================
    % Picard options
    %==============================================================
    picardOptions = make_Picard_options( ...
        input.Picard, verbose);

    %==============================================================
    % q-evaluation options
    %==============================================================
    [qInteriorPsiN, qContourNR] = ...
        make_q_evaluation_options(input);

    %==============================================================
    % A.1 Generate candidate interior P1 points
    %==============================================================
    meshTimer = tic;

    Rgrid = linspace( ...
        min(R_bnd), max(R_bnd), NR);

    Zgrid = linspace( ...
        min(Z_bnd), max(Z_bnd), NZ);

    [RR, ZZ] = meshgrid(Rgrid, Zgrid);

    [insideDomain, onBoundary] = inpolygon( ...
        RR, ZZ, R_bnd, Z_bnd);

    strictlyInside = ...
        insideDomain & ~onBoundary;

    interiorPoints = [ ...
        RR(strictlyInside), ...
        ZZ(strictlyInside)];

    %==============================================================
    % Include the nominal geometric center
    %==============================================================
    centerPoint = [R0, 0.0];

    [centerInside, centerOnBoundary] = inpolygon( ...
        centerPoint(1), centerPoint(2), ...
        R_bnd, Z_bnd);

    if ~(centerInside || centerOnBoundary)
        error('GS:PrepareProblem:CenterOutsideBoundary', ...
            ['The nominal center point [R0,0] lies outside ', ...
             'the prescribed plasma boundary.']);
    end

    distanceFromCenter = hypot( ...
        interiorPoints(:,1)-centerPoint(1), ...
        interiorPoints(:,2)-centerPoint(2));

    centerTolerance = ...
        1.e-12*coordinateScale;

    interiorPoints( ...
        distanceFromCenter < centerTolerance, :) = [];

    plasmaBoundary = [R_bnd, Z_bnd];

    P1points = [ ...
        plasmaBoundary; ...
        interiorPoints; ...
        centerPoint];

    boundaryConstraints = [ ...
        (1:nBoundary).', ...
        [2:nBoundary, 1].'];

    %==============================================================
    % A.2 Generate constrained P1 mesh
    %==============================================================
    meshP1 = generate_P1_mesh( ...
        P1points, boundaryConstraints, false);

    triangleCenters = incenter(meshP1);

    triangleInside = inpolygon( ...
        triangleCenters(:,1), ...
        triangleCenters(:,2), ...
        R_bnd, Z_bnd);

    if any(~triangleInside)
        error('GS:PrepareProblem:TriangleOutsideBoundary', ...
            ['The constrained triangulation contains triangles ', ...
             'outside the prescribed plasma boundary.']);
    end

    P1meshTime = toc(meshTimer);

    %==============================================================
    % A.3 Upgrade to P3 solution mesh
    %==============================================================
    P3timer = tic;

    P3 = generate_P3_nodes(meshP1);

    P3meshTime = toc(P3timer);

    Ndof = size(P3.points, 1);
    Nt = size(P3.elements, 1);

    if any(P3.points(:,1) <= 0)
        error('GS:PrepareProblem:NonpositiveRadius', ...
            'The P3 mesh contains R <= 0.');
    end

    if isempty(P3.boundaryNodes)
        error('GS:PrepareProblem:BoundaryNodes', ...
            'The P3 mesh contains no boundary nodes.');
    end

    % A reusable element locator for post-processing.
    P1locator = triangulation( ...
        P3.P1elements, P3.P1points);

    %==============================================================
    % A.4 Precompute P3 quadrature
    %==============================================================
    quadratureTimer = tic;

    quad = precompute_quadrature( ...
        quadratureDegree);

    quadratureTime = toc(quadratureTimer);

    %==============================================================
    % A.5 Assemble fixed stiffness matrix
    %==============================================================
    operatorTimer = tic;

    K = assemble_K_matrix(P3, quad);

    if ~issparse(K)
        warning('GS:PrepareProblem:NonSparseK', ...
            'The global stiffness matrix is not sparse.');
    end

    if any(~isfinite(nonzeros(K)))
        error('GS:PrepareProblem:InvalidK', ...
            'The global stiffness matrix contains NaN or Inf.');
    end

    %==============================================================
    % A.6 Apply fixed Dirichlet partition
    %==============================================================
    zeroRHS = zeros(Ndof, 1);

    [KII, ~, bc] = apply_BC( ...
        P3, K, zeroRHS, psiBoundary);

    freeNodes = bc.freeNodes;
    boundaryNodes = bc.boundaryNodes;
    boundaryValues = bc.boundaryValues;

    KIB = K(freeNodes, boundaryNodes);

    boundaryCorrection = ...
        KIB*boundaryValues;

    symmetryError = ...
        norm(KII-KII.', 'fro') ...
        /max(norm(KII, 'fro'), eps);

    if symmetryError > ...
            picardOptions.symmetryTolerance
        error('GS:PrepareProblem:Symmetry', ...
            ['KII symmetry error %.3e exceeds ', ...
             'the tolerance %.3e.'], ...
            symmetryError, ...
            picardOptions.symmetryTolerance);
    end

    %==============================================================
    % A.7 Factorize the fixed reduced matrix
    %==============================================================
    try
        Kfactor = decomposition( ...
            KII, 'chol', 'lower');
    catch ME
        error('GS:PrepareProblem:Factorization', ...
            ['Failed to construct the Cholesky ', ...
             'factorization of KII:\n%s'], ...
            ME.message);
    end

    operatorTime = toc(operatorTimer);

    %==============================================================
    % B.1 Construct or load the initial psi
    %==============================================================
    [psi0, initialInfo] = construct_initial_psi( ...
        input, P3, R_bnd, Z_bnd, ...
        psiBoundary, coordinateScale);

    %==============================================================
    % Common normalized-profile grid
    %==============================================================
    psiNProfile = ...
        linspace(0, 1, nProfile).';

    %==============================================================
    % Package prepared problem
    %==============================================================
    problem = struct();

    problem.prepared = true;

    % Mesh and interpolation geometry
    problem.meshP1 = meshP1;
    problem.mesh = P3;
    problem.P1locator = P1locator;

    % FEM quadrature
    problem.quad = quad;

    % Fixed linear operator
    problem.operator = struct();

    problem.operator.K = K;
    problem.operator.KII = KII;
    problem.operator.KIB = KIB;
    problem.operator.Kfactor = Kfactor;

    problem.operator.bc = bc;

    problem.operator.freeNodes = freeNodes;
    problem.operator.boundaryNodes = boundaryNodes;
    problem.operator.boundaryValues = boundaryValues;
    problem.operator.boundaryCorrection = ...
        boundaryCorrection;

    problem.operator.symmetryError = ...
        symmetryError;

    % Boundary information
    problem.boundary = struct();

    problem.boundary.R = R_bnd;
    problem.boundary.Z = Z_bnd;
    problem.boundary.points = plasmaBoundary;
    problem.boundary.psi = psiBoundary;
    problem.boundary.nPoint = nBoundary;

    % Dimensional parameters
    problem.dim = struct();

    problem.dim.R0 = R0;
    problem.dim.B0 = B0;
    problem.dim.aMinor = aMinor;

    % Grid information
    problem.grid = struct();

    problem.grid.NR = NR;
    problem.grid.NZ = NZ;
    problem.grid.nBoundary = nBoundary;
    problem.grid.nProfile = nProfile;
    problem.grid.quadratureDegree = ...
        quadratureDegree;

    problem.grid.Ndof = Ndof;
    problem.grid.Nelement = Nt;

    % Profile grid
    problem.profile = struct();
    problem.profile.psiN = psiNProfile;

    % Initial state
    problem.initial = initialInfo;
    problem.initial.psi = psi0;

    % Solver options
    problem.options = struct();
    problem.options.Picard = picardOptions;

    % q-evaluation settings
    problem.qEvaluation = struct();

    problem.qEvaluation.interiorPsiN = ...
        qInteriorPsiN;

    problem.qEvaluation.contourNR = ...
        qContourNR;

    % Output options
    problem.output = struct();

    problem.output.verbose = verbose;
    problem.output.checkTime = checkTime;
    problem.output.showPlot = showPlot;

    % Timings
    problem.timing = struct();

    problem.timing.P1mesh = P1meshTime;
    problem.timing.P3mesh = P3meshTime;
    problem.timing.quadrature = quadratureTime;
    problem.timing.operator = operatorTime;

    %==============================================================
    % Report preparation
    %==============================================================
    if verbose
        fprintf('\n');
        fprintf('=======================================================\n');
        fprintf(' Prepared Grad-Shafranov problem\n');
        fprintf('=======================================================\n');

        fprintf('P1 vertices          = %d\n', ...
            size(P3.P1points,1));

        fprintf('P1 triangles         = %d\n', Nt);
        fprintf('P3 DOFs              = %d\n', Ndof);

        fprintf('P3 boundary DOFs     = %d\n', ...
            numel(boundaryNodes));

        fprintf('P3 free DOFs         = %d\n', ...
            numel(freeNodes));

        fprintf('Quadrature degree    = %d\n', ...
            quadratureDegree);

        fprintf('Quadrature points    = %d\n', ...
            numel(quad.w));

        fprintf('KII symmetry error   = %.3e\n', ...
            symmetryError);

        fprintf('=======================================================\n\n');
    end

    if checkTime
        fprintf('P1 mesh time         = %.3f s\n', ...
            P1meshTime);

        fprintf('P3 mesh time         = %.3f s\n', ...
            P3meshTime);

        fprintf('Quadrature time      = %.3f s\n', ...
            quadratureTime);

        fprintf('Operator setup time  = %.3f s\n', ...
            operatorTime);
    end
end


function options = make_Picard_options( ...
    inputOptions, verbose)
%MAKE_PICARD_OPTIONS Convert input.Picard into the existing option format.

    options = struct();

    options.omega = ...
        inputOptions.omega;

    options.maxIterations = ...
        inputOptions.maxIterations;

    options.updateTolerance = ...
        inputOptions.updateTolerance;

    options.residualTolerance = ...
        inputOptions.residualTolerance;

    options.axisMode = ...
        inputOptions.axisMode;

    options.symmetryTolerance = ...
        inputOptions.symmetryTolerance;

    options.fluxSpanTolerance = ...
        inputOptions.fluxSpanTolerance;

    options.verbose = verbose;

    validateattributes(options.omega, {'numeric'}, ...
        {'real', 'finite', 'scalar', '>', 0, '<=', 1});

    validateattributes(options.maxIterations, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', 'positive'});

    validateattributes(options.updateTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(options.residualTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(options.symmetryTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    validateattributes(options.fluxSpanTolerance, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'});

    if isstring(options.axisMode)
        if ~isscalar(options.axisMode)
            error('GS:PrepareProblem:AxisMode', ...
                'axisMode must be scalar.');
        end

        options.axisMode = char(options.axisMode);
    end

    if ~ischar(options.axisMode) || ...
            ~isrow(options.axisMode)
        error('GS:PrepareProblem:AxisMode', ...
            'axisMode must be ''max'', ''min'', or ''auto''.');
    end

    options.axisMode = ...
        lower(strtrim(options.axisMode));

    if ~ismember(options.axisMode, ...
            {'max', 'min', 'auto'})
        error('GS:PrepareProblem:AxisMode', ...
            'axisMode must be ''max'', ''min'', or ''auto''.');
    end
end


function [interiorPsiN, contourNR] = ...
    make_q_evaluation_options(input)
%MAKE_Q_EVALUATION_OPTIONS Parse optional flux-surface sampling controls.

    interiorPsiN = unique([ ...
        0.05:0.05:0.80, ...
        0.825:0.025:0.95]).';
    contourNR = 500;

    if ~isfield(input, 'qEvaluation') || ...
            isempty(input.qEvaluation)
        return
    end

    options = input.qEvaluation;

    if ~isstruct(options) || ~isscalar(options)
        error('GS:PrepareProblem:QEvaluationType', ...
            'input.qEvaluation must be a scalar structure.');
    end

    allowedFields = {'interiorPsiN', 'contourNR'};
    unknownFields = setdiff(fieldnames(options), allowedFields);

    if ~isempty(unknownFields)
        error('GS:PrepareProblem:QEvaluationOption', ...
            'Unknown input.qEvaluation option: %s', ...
            strjoin(unknownFields, ', '));
    end

    if isfield(options, 'interiorPsiN') && ...
            ~isempty(options.interiorPsiN)
        interiorPsiN = options.interiorPsiN(:);
    end

    if numel(interiorPsiN) < 3 || ...
            any(~isfinite(interiorPsiN)) || ...
            any(interiorPsiN <= 0) || ...
            any(interiorPsiN >= 1) || ...
            any(diff(interiorPsiN) <= 0)
        error('GS:PrepareProblem:QEvaluationLevels', ...
            ['input.qEvaluation.interiorPsiN must contain at least ', ...
             'three finite, strictly increasing values in (0,1).']);
    end

    if isfield(options, 'contourNR') && ...
            ~isempty(options.contourNR)
        contourNR = options.contourNR;
    end

    validateattributes(contourNR, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2}, ...
        mfilename, 'input.qEvaluation.contourNR');
end


function [psi0, info] = construct_initial_psi( ...
    input, P3, R_bnd, Z_bnd, ...
    psiBoundary, coordinateScale)
%CONSTRUCT_INITIAL_PSI
%
% Load an explicitly supplied psi initial condition or reproduce the
% distance-to-boundary initial condition currently used by solve_GS.

    Ndof = size(P3.points,1);

    hasInitialPsi = ...
        isfield(input, 'initial') && ...
        isstruct(input.initial) && ...
        isfield(input.initial, 'psi') && ...
        ~isempty(input.initial.psi);

    if hasInitialPsi
        psi0 = input.initial.psi(:);

        if numel(psi0) ~= Ndof
            error('GS:PrepareProblem:InitialPsiSize', ...
                ['input.initial.psi must contain one value ', ...
                 'per P3 degree of freedom.']);
        end

        if any(~isfinite(psi0))
            error('GS:PrepareProblem:InitialPsi', ...
                'input.initial.psi contains NaN or Inf.');
        end

        initialBoundaryError = norm( ...
            psi0(P3.boundaryNodes)-psiBoundary, inf);

        psi0(P3.boundaryNodes) = psiBoundary;

        info = struct();

        info.type = 'user-supplied';
        info.psiAxisGuess = NaN;
        info.initialBoundaryError = ...
            initialBoundaryError;

        return
    end

    % Preserve the current solve_GS default.
    psiAxisGuess = 0.5;

    if isfield(input, 'initial') && ...
            isstruct(input.initial) && ...
            isfield(input.initial, 'psiAxisGuess') && ...
            ~isempty(input.initial.psiAxisGuess)

        psiAxisGuess = ...
            input.initial.psiAxisGuess;
    end

    validateattributes(psiAxisGuess, {'numeric'}, ...
        {'real', 'finite', 'scalar'});

    fluxGuessScale = max([ ...
        abs(psiBoundary), ...
        abs(psiAxisGuess), ...
        1]);

    if abs(psiAxisGuess-psiBoundary) <= ...
            100*eps(fluxGuessScale)
        error('GS:PrepareProblem:InitialFluxSpan', ...
            ['psiAxisGuess must differ from the prescribed ', ...
             'boundary flux.']);
    end

    distanceToBoundary = inf(Ndof,1);

    nBoundary = numel(R_bnd);

    for k = 1:nBoundary
        distanceToPoint = hypot( ...
            P3.points(:,1)-R_bnd(k), ...
            P3.points(:,2)-Z_bnd(k));

        distanceToBoundary = min( ...
            distanceToBoundary, distanceToPoint);
    end

    maximumDistance = ...
        max(distanceToBoundary);

    if ~isfinite(maximumDistance) || ...
            maximumDistance <= ...
            100*eps(coordinateScale)
        error('GS:PrepareProblem:InitialPsi', ...
            'Failed to construct the initial flux guess.');
    end

    normalizedDistance = ...
        distanceToBoundary/maximumDistance;

    psi0 = psiBoundary ...
        +(psiAxisGuess-psiBoundary)*normalizedDistance;

    initialBoundaryError = norm( ...
        psi0(P3.boundaryNodes)-psiBoundary, inf);

    psi0(P3.boundaryNodes) = psiBoundary;

    info = struct();

    info.type = 'distance-to-boundary';
    info.psiAxisGuess = psiAxisGuess;
    info.initialBoundaryError = ...
        initialBoundaryError;
end


function assert_required_fields( ...
    value, requiredFields, valueName)
%ASSERT_REQUIRED_FIELDS Check required structure fields.

    if ~isstruct(value) || ~isscalar(value)
        error('GS:PrepareProblem:StructureType', ...
            '%s must be a scalar structure.', valueName);
    end

    for k = 1:numel(requiredFields)
        fieldName = requiredFields{k};

        if ~isfield(value, fieldName)
            error('GS:PrepareProblem:MissingField', ...
                '%s.%s is required.', ...
                valueName, fieldName);
        end
    end
end
