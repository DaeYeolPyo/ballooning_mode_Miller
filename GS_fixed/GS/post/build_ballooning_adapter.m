function [eq, eqfunc, diagnostics] = build_ballooning_adapter( ...
    equilibrium, opts)
%BUILD_BALLOONING_ADAPTER Expose a P3 GS equilibrium to legacy consumers.
%
%   [eq, eqfunc, diagnostics] = BUILD_BALLOONING_ADAPTER(equilibrium)
%
% The returned objects follow the interface produced by read_geqdsk and
% build_interpolants.  A Cartesian psi grid is generated only for contour
% discovery.  eqfunc.psi, eqfunc.psiR, and eqfunc.psiZ evaluate the P3 FEM
% field directly and therefore do not differentiate the Cartesian cache.
%
% The adapter uses an outward-increasing flux coordinate
%
%   psi_adapter = abs(nativeSpan)*psiN,
%
% so its magnetic axis is psi=0 and its fixed boundary is
% psi=abs(nativeSpan).  Profile derivatives are converted from d/dpsiN to
% d/dpsi_adapter.  The native GS convention is retained in diagnostics.
%
% Name-value options
%   NR : number of Cartesian R points and legacy profile samples (257)
%   NZ : number of Cartesian Z points.  NaN selects an aspect-ratio value.

    arguments
        equilibrium (1,1) struct
        opts.NR (1,1) double = 257
        opts.NZ (1,1) double = NaN
    end

    validateattributes(opts.NR, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 16}, ...
        mfilename, 'NR');

    requiredEquilibriumFields = { ...
        'mesh', 'psi', 'profile', 'global', 'qResult'};
    assert_fields(equilibrium, requiredEquilibriumFields, 'equilibrium');

    requiredProfileFields = { ...
        'psiN', 'pressure', 'dp_dpsiN', ...
        'Fpol', 'FdF_dpsiN', 'q'};
    assert_fields(equilibrium.profile, requiredProfileFields, ...
        'equilibrium.profile');

    requiredGlobalFields = {'simag', 'sibry', 'Ip'};
    assert_fields(equilibrium.global, requiredGlobalFields, ...
        'equilibrium.global');

    P3 = equilibrium.mesh;
    requiredMeshFields = { ...
        'points', 'elements', 'P1points', 'P1elements'};
    assert_fields(P3, requiredMeshFields, 'equilibrium.mesh');

    nativePsi = equilibrium.psi(:);
    if numel(nativePsi) ~= size(P3.points, 1) || ...
            any(~isfinite(nativePsi))
        error('GS:BallooningAdapter:InvalidPsi', ...
            ['equilibrium.psi must contain one finite value per ', ...
             'global P3 degree of freedom.']);
    end

    nativeAxis = equilibrium.global.simag;
    nativeBoundary = equilibrium.global.sibry;
    nativeSpan = nativeBoundary-nativeAxis;

    if ~isfinite(nativeSpan) || nativeSpan == 0
        error('GS:BallooningAdapter:InvalidFluxSpan', ...
            'The native boundary-axis flux span must be finite and nonzero.');
    end

    outwardSpan = abs(nativeSpan);
    spatialScale = outwardSpan/nativeSpan;

    profilePsiN = equilibrium.profile.psiN(:);
    validate_profile_grid(profilePsiN);

    profileNames = { ...
        'pressure', 'dp_dpsiN', 'Fpol', 'FdF_dpsiN', 'q'};
    for k = 1:numel(profileNames)
        name = profileNames{k};
        value = equilibrium.profile.(name);
        if numel(value) ~= numel(profilePsiN) || any(~isfinite(value))
            error('GS:BallooningAdapter:InvalidProfile', ...
                ['equilibrium.profile.%s must contain one finite ', ...
                 'value per profile psiN sample.'], name);
        end
    end

    boundary = get_ordered_boundary(equilibrium, P3);
    [R0, B0, aMinor, usedReferenceFallback] = ...
        get_reference_values(equilibrium, boundary);

    axisPoint = get_axis_point(equilibrium);

    Rlimits = [min(P3.P1points(:,1)), max(P3.P1points(:,1))];
    Zlimits = [min(P3.P1points(:,2)), max(P3.P1points(:,2))];

    Rspan = diff(Rlimits);
    Zspan = diff(Zlimits);
    if Rspan <= 0 || Zspan <= 0
        error('GS:BallooningAdapter:InvalidDomain', ...
            'The FEM domain must have positive R and Z extents.');
    end

    nR = double(opts.NR);
    if isnan(opts.NZ)
        nZ = max(16, round(nR*Zspan/Rspan));
    else
        validateattributes(opts.NZ, {'numeric'}, ...
            {'real', 'finite', 'scalar', 'integer', '>=', 16}, ...
            mfilename, 'NZ');
        nZ = double(opts.NZ);
    end

    rgrid = linspace(Rlimits(1), Rlimits(2), nR);
    zgrid = linspace(Zlimits(1), Zlimits(2), nZ);
    [Rcache, Zcache] = ndgrid(rgrid, zgrid);

    P1 = triangulation(P3.P1elements, P3.P1points);

    valueCache = build_P3_value_cache(P3,[Rcache(:), Zcache(:)]);
    nativeGrid = valueCache.operator*nativePsi;
    nativeGrid(~valueCache.inside)=NaN;
    psiNGrid = (nativeGrid-nativeAxis)/nativeSpan;
    adapterGrid = outwardSpan*psiNGrid;
    psirz = reshape(adapterGrid, [nR, nZ]);

    % GEQDSK uses nw both as the R-grid size and the profile-array length.
    % Resample all profiles to nR so legacy fit functions keep that contract.
    adapterPsiN = linspace(0, 1, nR).';

    pressure = interp1(profilePsiN, ...
        equilibrium.profile.pressure(:), adapterPsiN, 'pchip');
    dpDpsiN = interp1(profilePsiN, ...
        equilibrium.profile.dp_dpsiN(:), adapterPsiN, 'pchip');
    Fpol = interp1(profilePsiN, ...
        equilibrium.profile.Fpol(:), adapterPsiN, 'pchip');
    FdFDpsiN = interp1(profilePsiN, ...
        equilibrium.profile.FdF_dpsiN(:), adapterPsiN, 'pchip');
    qpsi = interp1(profilePsiN, ...
        equilibrium.profile.q(:), adapterPsiN, 'pchip');

    pprime = dpDpsiN/outwardSpan;
    ffprim = FdFDpsiN/outwardSpan;

    eq = struct();
    eq.filename = '';
    eq.header = 'P3 FEM GS ballooning adapter';
    eq.idum = 0;

    eq.nw = nR;
    eq.nh = nZ;

    eq.rdim = Rspan;
    eq.zdim = Zspan;
    eq.rcentr = R0;
    eq.rleft = Rlimits(1);
    eq.zmid = mean(Zlimits);

    eq.rmaxis = axisPoint(1);
    eq.zmaxis = axisPoint(2);

    eq.simag = 0;
    eq.sibry = outwardSpan;
    eq.psi_axis = eq.simag;
    eq.psi_bdry = eq.sibry;

    eq.bcentr = B0;
    eq.current = equilibrium.global.Ip;

    eq.fpol = Fpol;
    eq.pres = pressure;
    eq.ffprim = ffprim;
    eq.pprime = pprime;
    eq.qpsi = qpsi;

    eq.rgrid = rgrid;
    eq.zgrid = zgrid;
    eq.psirz = psirz;
    eq.psin = adapterPsiN.';

    eq.nbbbs = size(boundary, 1);
    eq.rbbbs = boundary(:, 1);
    eq.zbbbs = boundary(:, 2);

    eq.limitr = 0;
    eq.rlim = [];
    eq.zlim = [];

    FInterpolator = griddedInterpolant( ...
        adapterPsiN, Fpol, 'linear', 'none');
    pressureInterpolator = griddedInterpolant( ...
        adapterPsiN, pressure, 'linear', 'none');
    pprimeInterpolator = griddedInterpolant( ...
        adapterPsiN, pprime, 'linear', 'none');
    FFprimeInterpolator = griddedInterpolant( ...
        adapterPsiN, ffprim, 'linear', 'none');

    qPiecewisePolynomial = pchip(adapterPsiN, qpsi);
    qprimeNPiecewisePolynomial = ...
        differentiate_scalar_pp(qPiecewisePolynomial);

    eqfunc = struct();

    eqfunc.psi = @(Rq, Zq) evaluate_spatial_component( ...
        P1, P3, nativePsi, nativeAxis, nativeSpan, ...
        outwardSpan, Rq, Zq, 'psi');

    eqfunc.psiR = @(Rq, Zq) evaluate_spatial_component( ...
        P1, P3, nativePsi, nativeAxis, nativeSpan, ...
        outwardSpan, Rq, Zq, 'psiR');

    eqfunc.psiZ = @(Rq, Zq) evaluate_spatial_component( ...
        P1, P3, nativePsi, nativeAxis, nativeSpan, ...
        outwardSpan, Rq, Zq, 'psiZ');

    eqfunc.psiN_from_psi = @(psi) psi/outwardSpan;
    eqfunc.psi_from_psiN = @(psiN) outwardSpan*psiN;

    eqfunc.F = FInterpolator;
    eqfunc.p = pressureInterpolator;
    eqfunc.pprime = pprimeInterpolator;
    eqfunc.FFprime = FFprimeInterpolator;
    eqfunc.q = @(psiN) ppval(qPiecewisePolynomial, psiN);
    eqfunc.qprimeN = @(psiN) ...
        ppval(qprimeNPiecewisePolynomial, psiN);
    eqfunc.qprime = @(psiN) ...
        ppval(qprimeNPiecewisePolynomial, psiN)/outwardSpan;

    eqfunc.eval = @(Rq, Zq) evaluate_equilibrium( ...
        eqfunc, P1, P3, nativePsi, nativeAxis, nativeSpan, ...
        outwardSpan, Rq, Zq);

    diagnostics = struct();
    diagnostics.nativeAxis = nativeAxis;
    diagnostics.nativeBoundary = nativeBoundary;
    diagnostics.nativeSpan = nativeSpan;
    diagnostics.outwardSpan = outwardSpan;
    diagnostics.spatialScale = spatialScale;
    diagnostics.gridSize = [nR, nZ];
    diagnostics.finiteGridFraction = nnz(isfinite(psirz))/numel(psirz);
    diagnostics.usedReferenceFallback = usedReferenceFallback;
    diagnostics.reference = struct( ...
        'R0', R0, 'B0', B0, 'aMinor', aMinor);
    diagnostics.fluxConvention = ...
        'psi_adapter=abs(nativeSpan)*psiN';

    eq.adapter = diagnostics;
end


function value = evaluate_spatial_component( ...
    P1, P3, nativeField, nativeAxis, nativeSpan, outwardSpan, ...
    Rq, Zq, component)

    [psi, psiR, psiZ] = evaluate_spatial_field( ...
        P1, P3, nativeField, nativeAxis, nativeSpan, ...
        outwardSpan, Rq, Zq);

    switch component
        case 'psi'
            value = psi;
        case 'psiR'
            value = psiR;
        case 'psiZ'
            value = psiZ;
        otherwise
            error('GS:BallooningAdapter:InternalComponent', ...
                'Unknown spatial component %s.', component);
    end
end


function [psi, psiR, psiZ] = evaluate_spatial_field( ...
    P1, P3, nativeField, nativeAxis, nativeSpan, outwardSpan, Rq, Zq)

    if ~isnumeric(Rq) || ~isnumeric(Zq) || ~isequal(size(Rq), size(Zq))
        error('GS:BallooningAdapter:QuerySize', ...
            'Rq and Zq must be numeric arrays with equal sizes.');
    end

    originalSize = size(Rq);
    query = [Rq(:), Zq(:)];

    [nativeValue, nativeGradient] = ...
        eval_P3_field(P1, P3, nativeField, query);

    spatialScale = outwardSpan/nativeSpan;

    adapterValue = outwardSpan*(nativeValue-nativeAxis)/nativeSpan;
    adapterGradient = spatialScale*nativeGradient;

    psi = reshape(adapterValue, originalSize);
    psiR = reshape(adapterGradient(:, 1), originalSize);
    psiZ = reshape(adapterGradient(:, 2), originalSize);
end


function out = evaluate_equilibrium( ...
    eqfunc, P1, P3, nativeField, nativeAxis, nativeSpan, ...
    outwardSpan, Rq, Zq)

    [psi, psiR, psiZ] = evaluate_spatial_field( ...
        P1, P3, nativeField, nativeAxis, nativeSpan, ...
        outwardSpan, Rq, Zq);

    out = struct();

    out.R = Rq;
    out.Z = Zq;
    out.psi = psi;
    out.psiR = psiR;
    out.psiZ = psiZ;
    out.psiN = psi/outwardSpan;

    out.gradPsi = hypot(psiR, psiZ);

    out.F = eqfunc.F(out.psiN);
    out.q = eqfunc.q(out.psiN);
    out.p = eqfunc.p(out.psiN);
    out.pprime = eqfunc.pprime(out.psiN);
    out.ffprime = eqfunc.FFprime(out.psiN);
    out.qprime = eqfunc.qprime(out.psiN);

    out.Bp = out.gradPsi./Rq;
    out.Bphi = out.F./Rq;
    out.B2 = out.Bp.^2 + out.Bphi.^2;
end


function point = get_axis_point(equilibrium)
    point = [];

    if isfield(equilibrium, 'qResult') && ...
            isstruct(equilibrium.qResult) && ...
            isfield(equilibrium.qResult, 'axis') && ...
            isstruct(equilibrium.qResult.axis) && ...
            isfield(equilibrium.qResult.axis, 'point')
        point = equilibrium.qResult.axis.point(:).';
    end

    if isempty(point) && isfield(equilibrium, 'axisPoint')
        point = equilibrium.axisPoint(:).';
    end

    if numel(point) ~= 2 || any(~isfinite(point)) || point(1) <= 0
        error('GS:BallooningAdapter:InvalidAxis', ...
            'A finite positive-R magnetic-axis point is required.');
    end
end


function boundary = get_ordered_boundary(equilibrium, P3)
    boundary = [];

    if isfield(equilibrium, 'boundary') && ...
            isstruct(equilibrium.boundary)
        if isfield(equilibrium.boundary, 'points') && ...
                ~isempty(equilibrium.boundary.points)
            boundary = equilibrium.boundary.points;
        elseif isfield(equilibrium.boundary, 'R') && ...
                isfield(equilibrium.boundary, 'Z')
            boundary = [ ...
                equilibrium.boundary.R(:), ...
                equilibrium.boundary.Z(:)];
        end
    end

    if isempty(boundary)
        if ~isfield(P3, 'boundaryEdges') || isempty(P3.boundaryEdges)
            error('GS:BallooningAdapter:MissingBoundary', ...
                ['The equilibrium has no ordered boundary and the ', ...
                 'mesh has no boundaryEdges fallback.']);
        end
        vertexOrder = order_boundary_edges(P3.boundaryEdges);
        boundary = P3.P1points(vertexOrder, :);
    end

    if size(boundary, 2) ~= 2 || size(boundary, 1) < 3 || ...
            any(~isfinite(boundary), 'all')
        error('GS:BallooningAdapter:InvalidBoundary', ...
            'The fixed boundary must be a finite N-by-2 array.');
    end

    scale = max([range(boundary(:, 1)), range(boundary(:, 2)), 1]);
    tolerance = 1.e-12*scale;

    keep = [true; hypot(diff(boundary(:, 1)), ...
        diff(boundary(:, 2))) > tolerance];
    boundary = boundary(keep, :);

    if size(boundary, 1) > 1 && ...
            norm(boundary(end, :)-boundary(1, :)) <= tolerance
        boundary(end, :) = [];
    end

    signedArea = 0.5*sum( ...
        boundary(:, 1).*circshift(boundary(:, 2), -1) ...
        - circshift(boundary(:, 1), -1).*boundary(:, 2));

    if signedArea < 0
        boundary = flipud(boundary);
    end
end


function order = order_boundary_edges(edges)
    edges = double(edges);
    vertices = unique(edges(:));

    degree = accumarray(edges(:), 1, [max(vertices), 1]);
    if any(degree(vertices) ~= 2)
        error('GS:BallooningAdapter:BoundaryTopology', ...
            'The mesh boundary must be one closed non-branching loop.');
    end

    used = false(size(edges, 1), 1);
    order = zeros(size(edges, 1), 1);

    order(1) = edges(1, 1);
    current = order(1);

    for k = 1:size(edges, 1)
        candidate = find(~used & ...
            (edges(:, 1) == current | edges(:, 2) == current), 1);

        if isempty(candidate)
            error('GS:BallooningAdapter:BoundaryTopology', ...
                'Failed to traverse the complete mesh boundary.');
        end

        used(candidate) = true;
        edge = edges(candidate, :);
        next = edge(1);
        if next == current
            next = edge(2);
        end

        if k < size(edges, 1)
            order(k+1) = next;
        elseif next ~= order(1)
            error('GS:BallooningAdapter:BoundaryTopology', ...
                'The traversed mesh boundary is not closed.');
        end

        current = next;
    end
end


function [R0, B0, aMinor, usedFallback] = ...
    get_reference_values(equilibrium, boundary)

    usedFallback = false;

    if isfield(equilibrium, 'reference') && ...
            isstruct(equilibrium.reference) && ...
            all(isfield(equilibrium.reference, ...
                {'R0', 'B0', 'aMinor'}))
        R0 = equilibrium.reference.R0;
        B0 = equilibrium.reference.B0;
        aMinor = equilibrium.reference.aMinor;
    else
        usedFallback = true;

        axisPoint = get_axis_point(equilibrium);
        R0 = axisPoint(1);
        aMinor = 0.5*range(boundary(:, 1));
        B0 = equilibrium.profile.Fpol(end)/R0;
    end

    if ~isfinite(R0) || R0 <= 0 || ...
            ~isfinite(B0) || B0 == 0 || ...
            ~isfinite(aMinor) || aMinor <= 0
        error('GS:BallooningAdapter:InvalidReference', ...
            'R0, B0, and aMinor reference values are invalid.');
    end
end


function validate_profile_grid(psiN)
    tolerance = 100*eps(max(1, max(abs(psiN))));

    if numel(psiN) < 2 || any(~isfinite(psiN)) || ...
            any(diff(psiN) <= 0) || ...
            abs(psiN(1)) > tolerance || ...
            abs(psiN(end)-1) > tolerance
        error('GS:BallooningAdapter:InvalidProfileGrid', ...
            ['equilibrium.profile.psiN must be finite, strictly ', ...
             'increasing, and span [0,1].']);
    end
end


function dpp = differentiate_scalar_pp(pp)
    [breaks, coefficients, ~, order, dimension] = unmkpp(pp);

    if dimension ~= 1
        error('GS:BallooningAdapter:NonScalarQ', ...
            'The q profile must be scalar-valued.');
    end

    if order <= 1
        derivativeCoefficients = zeros(size(coefficients, 1), 1);
    else
        powers = order-1:-1:1;
        derivativeCoefficients = ...
            coefficients(:, 1:end-1).*powers;
    end

    dpp = mkpp(breaks, derivativeCoefficients);
end


function assert_fields(value, names, variableName)
    if ~isstruct(value) || ~isscalar(value)
        error('GS:BallooningAdapter:InvalidStructure', ...
            '%s must be a scalar structure.', variableName);
    end

    for k = 1:numel(names)
        name = names{k};
        if ~isfield(value, name)
            error('GS:BallooningAdapter:MissingField', ...
                '%s.%s is required.', variableName, name);
        end
    end
end
