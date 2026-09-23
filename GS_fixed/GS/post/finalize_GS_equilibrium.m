function equilibrium = finalize_GS_equilibrium(problem, solved)
%FINALIZE_GS_EQUILIBRIUM
% Construct the public equilibrium output from a completed GS solve.
%
%   equilibrium = finalize_GS_equilibrium(problem, solved)
%
% This function performs only mode-independent postprocessing:
%
%   1. pressure reconstruction from pprime
%   2. plasma-current evaluation
%   3. q interpolation onto the common profile grid
%   4. magnetic-axis toroidal-field evaluation
%   5. global-parameter evaluation
%   6. equilibrium structure construction
%
% It must not perform:
%
%   - Picard iteration
%   - p-q outer iteration
%   - G-profile updating
%   - q-geometry reconstruction
%
% Required solved fields
% ----------------------
%
%   solved.mode
%   solved.state
%   solved.innerResult
%   solved.pprime
%   solved.pedge
%   solved.FProfile
%   solved.qResult
%
% Optional solved fields
% ----------------------
%
%   solved.solveTime
%   solved.constraintResult
%   solved.FFprimeModel
%   solved.qTargetEvaluated
%   solved.qGeometry

    narginchk(2, 2);

    %==============================================================
    % Input validation
    %==============================================================
    assert_required_fields(problem, ...
        {'mesh', 'quad', 'boundary', 'dim', 'profile'}, ...
        'problem');

    assert_required_fields(problem.boundary, ...
        {'psi'}, 'problem.boundary');

    assert_required_fields(problem.dim, ...
        {'R0', 'B0', 'aMinor'}, 'problem.dim');

    assert_required_fields(problem.profile, ...
        {'psiN'}, 'problem.profile');

    assert_required_fields(solved, ...
        {'mode', 'state', 'innerResult', ...
         'pprime', 'pedge', 'FProfile', 'qResult'}, ...
        'solved');

    state = solved.state;
    innerResult = solved.innerResult;
    FProfile = solved.FProfile;
    qResult = solved.qResult;

    assert_required_fields(state, ...
        {'psi', 'psiN', 'psiAxis', 'dpsi', ...
         'axisNode', 'axisPoint'}, ...
        'solved.state');

    assert_required_fields(innerResult, ...
        {'finalRHS'}, ...
        'solved.innerResult');

    assert_required_fields(FProfile, ...
        {'psiN', 'G', 'F', 'FFprime'}, ...
        'solved.FProfile');

    assert_required_fields(qResult, ...
        {'psiN', 'q', 'axis'}, ...
        'solved.qResult');

    assert_required_fields(qResult.axis, ...
        {'point'}, ...
        'solved.qResult.axis');

    profileMode = solved.mode;

    validateattributes(profileMode, {'numeric'}, ...
        {'scalar', 'integer'});

    if ~ismember(profileMode, [1, 3])
        error('GS:Finalize:UnsupportedMode', ...
            'Only profile modes 1 and 3 are supported.');
    end

    %==============================================================
    % Prepared problem data
    %==============================================================
    P3 = problem.mesh;
    quad = problem.quad;

    R0 = problem.dim.R0;
    B0 = problem.dim.B0;
    aMinor = problem.dim.aMinor;

    psiBoundary = problem.boundary.psi;
    psiNProfile = problem.profile.psiN(:);

    %==============================================================
    % Solution state
    %==============================================================
    psi = state.psi(:);
    psiN = state.psiN(:);

    psiAxis = state.psiAxis;
    dpsi = state.dpsi;

    axisNode = state.axisNode;
    axisPoint = state.axisPoint(:).';

    nDof = size(P3.points, 1);

    if numel(psi) ~= nDof
        error('GS:Finalize:PsiSize', ...
            'state.psi must contain one value per P3 node.');
    end

    if numel(psiN) ~= nDof
        error('GS:Finalize:PsiNSize', ...
            'state.psiN must contain one value per P3 node.');
    end

    if any(~isfinite(psi)) || any(~isfinite(psiN))
        error('GS:Finalize:InvalidFlux', ...
            'state.psi or state.psiN contains NaN or Inf.');
    end

    if numel(axisPoint) ~= 2 || any(~isfinite(axisPoint))
        error('GS:Finalize:InvalidAxisPoint', ...
            'state.axisPoint must contain finite [R,Z] coordinates.');
    end

    if ~isscalar(axisNode) || ...
            axisNode < 1 || axisNode > nDof
        error('GS:Finalize:InvalidAxisNode', ...
            'state.axisNode is invalid.');
    end

    %==============================================================
    % Common profile data
    %==============================================================
    pprime = solved.pprime(:);
    pedge = solved.pedge;

    if numel(pprime) ~= numel(psiNProfile)
        error('GS:Finalize:PprimeSize', ...
            ['solved.pprime must have the same length as ', ...
             'problem.profile.psiN.']);
    end

    if any(~isfinite(pprime))
        error('GS:Finalize:InvalidPprime', ...
            'solved.pprime contains NaN or Inf.');
    end

    validateattributes(pedge, {'numeric'}, ...
        {'real', 'finite', 'scalar'});

    profilePsiN = FProfile.psiN(:);
    Gpol = FProfile.G(:);
    Fpol = FProfile.F(:);
    FFprime = FProfile.FFprime(:);

    nProfile = numel(psiNProfile);

    if numel(profilePsiN) ~= nProfile || ...
            numel(Gpol) ~= nProfile || ...
            numel(Fpol) ~= nProfile || ...
            numel(FFprime) ~= nProfile
        error('GS:Finalize:FProfileSize', ...
            ['FProfile.psiN, G, F, and FFprime must use the ', ...
             'common problem.profile.psiN grid.']);
    end

    profileGridTolerance = ...
        100*eps(max(1, max(abs(psiNProfile))));

    if norm(profilePsiN - psiNProfile, inf) > ...
            profileGridTolerance
        error('GS:Finalize:FProfileGrid', ...
            ['solved.FProfile must be sampled on ', ...
             'problem.profile.psiN.']);
    end

    if any(~isfinite(Gpol)) || any(Gpol <= 0)
        error('GS:Finalize:InvalidG', ...
            'The final G profile must be finite and positive.');
    end

    if any(~isfinite(Fpol)) || any(~isfinite(FFprime))
        error('GS:Finalize:InvalidFProfile', ...
            'The final F or FFprime profile contains NaN or Inf.');
    end

    %==============================================================
    % Pressure reconstruction
    %==============================================================
    integratedP = cumtrapz(psiNProfile, pprime);

    pressure = ...
        pedge + integratedP - integratedP(end);

    if any(~isfinite(pressure))
        error('GS:Finalize:InvalidPressure', ...
            'The reconstructed pressure contains NaN or Inf.');
    end

    %==============================================================
    % Plasma current
    %==============================================================
    mu0 = 4*pi*1.e-7;

    finalRHS = innerResult.finalRHS(:);

    if numel(finalRHS) ~= nDof
        error('GS:Finalize:RHSSize', ...
            'innerResult.finalRHS must contain one value per P3 node.');
    end

    if any(~isfinite(finalRHS))
        error('GS:Finalize:InvalidRHS', ...
            'innerResult.finalRHS contains NaN or Inf.');
    end

    Ip = sum(finalRHS)/mu0;

    %==============================================================
    % q interpolation onto the common profile grid
    %==============================================================
    qPsiN = qResult.psiN(:);
    qValues = qResult.q(:);

    if numel(qPsiN) ~= numel(qValues)
        error('GS:Finalize:QSize', ...
            'qResult.psiN and qResult.q must have equal lengths.');
    end

    if numel(qPsiN) < 2 || ...
            any(~isfinite(qPsiN)) || ...
            any(diff(qPsiN) <= 0)
        error('GS:Finalize:QGrid', ...
            'qResult.psiN must be finite and strictly increasing.');
    end

    if any(~isfinite(qValues))
        error('GS:Finalize:InvalidQ', ...
            'qResult.q contains NaN or Inf.');
    end

    qGridTolerance = ...
        100*eps(max(1, max(abs(qPsiN))));

    if qPsiN(1) > qGridTolerance || ...
            qPsiN(end) < 1-qGridTolerance
        error('GS:Finalize:QRange', ...
            'qResult.psiN must cover the interval [0,1].');
    end

    qProfile = interp1( ...
        qPsiN, qValues, ...
        psiNProfile, 'pchip');

    if any(~isfinite(qProfile))
        error('GS:Finalize:QInterpolation', ...
            'Failed to interpolate q onto the common profile grid.');
    end

    %==============================================================
    % Magnetic-axis toroidal field
    %==============================================================
    qAxisPoint = qResult.axis.point(:).';

    if numel(qAxisPoint) ~= 2 || ...
            any(~isfinite(qAxisPoint)) || ...
            qAxisPoint(1) <= 0
        error('GS:Finalize:QAxisPoint', ...
            'qResult.axis.point must contain valid [R,Z] coordinates.');
    end

    Faxis = interp1( ...
        profilePsiN, Fpol, 0, 'linear');

    if ~isfinite(Faxis)
        error('GS:Finalize:FAxis', ...
            'Failed to evaluate F at the magnetic axis.');
    end

    Baxis = Faxis/qAxisPoint(1);

    %==============================================================
    % Global parameters
    %==============================================================
    globalResult = evaluate_global_parameters_local( ...
        P3, quad, psi, ...
        psiAxis, dpsi, ...
        psiNProfile, pressure, ...
        R0, B0, aMinor, Ip);

    %==============================================================
    % Public equilibrium structure
    %==============================================================
    equilibrium = struct();

    equilibrium.mesh = P3;

    equilibrium.psi = psi;
    equilibrium.psiN = psiN;

    equilibrium.axisNode = axisNode;
    equilibrium.axisPoint = axisPoint;

    % Profiles
    equilibrium.profile = struct();

    equilibrium.profile.mode = profileMode;
    equilibrium.profile.psiN = psiNProfile;

    equilibrium.profile.pressure = pressure;
    equilibrium.profile.dp_dpsiN = pprime;

    equilibrium.profile.Gpol = Gpol;
    equilibrium.profile.Fpol = Fpol;
    equilibrium.profile.FdF_dpsiN = FFprime;

    equilibrium.profile.q = qProfile;

    % Global parameters
    equilibrium.global = struct();

    equilibrium.global.Raxis = axisPoint(1);
    equilibrium.global.Zaxis = axisPoint(2);

    equilibrium.global.simag = psiAxis;
    equilibrium.global.sibry = psiBoundary;

    equilibrium.global.Baxis = Baxis;
    equilibrium.global.Ip = Ip;

    equilibrium.global.volume = ...
        globalResult.volume;

    equilibrium.global.betaTor = ...
        globalResult.betaTor;

    equilibrium.global.betaPol = ...
        globalResult.betaPol;

    equilibrium.global.li3 = ...
        globalResult.li3;

    equilibrium.global.betaN = ...
        globalResult.betaN;

    % Keep backward compatibility with report_GS_result and
    % plot_GS_equilibrium. For mode 3 this is the final inner Picard
    % result, not the p-q outer result.
    equilibrium.solverResult = innerResult;

    %==============================================================
    % Mode-3-specific results
    %==============================================================
    if profileMode == 3
        assert_required_fields(solved, ...
            {'constraintResult', 'FFprimeModel'}, ...
            'solved');

        equilibrium.constraintResult = ...
            solved.constraintResult;

        equilibrium.FFprimeModel = ...
            solved.FFprimeModel;

        if isfield(solved, 'qTargetEvaluated') && ...
                ~isempty(solved.qTargetEvaluated)

            qTargetEvaluated = ...
                solved.qTargetEvaluated(:);

            if numel(qTargetEvaluated) ~= numel(qPsiN)
                error('GS:Finalize:QTargetSize', ...
                    ['solved.qTargetEvaluated must have the same ', ...
                     'length as qResult.psiN.']);
            end

            equilibrium.profile.qPrescribed = struct();

            equilibrium.profile.qPrescribed.psiN = ...
                qPsiN;

            equilibrium.profile.qPrescribed.q = ...
                qTargetEvaluated;
        end
    end

    %==============================================================
    % Optional diagnostic data
    %==============================================================
    if isfield(solved, 'qGeometry') && ...
            ~isempty(solved.qGeometry)
        equilibrium.qGeometry = solved.qGeometry;
    end

    equilibrium.qResult = qResult;

    if isfield(problem, 'timing')
        equilibrium.timing = problem.timing;
    else
        equilibrium.timing = struct();
    end

    if isfield(solved, 'solveTime') && ...
            ~isempty(solved.solveTime)
        equilibrium.timing.solve = solved.solveTime;
    end
end


function result = evaluate_global_parameters_local( ...
    P3, quad, psi, ...
    psiAxis, dpsi, ...
    psiNTable, pressureTable, ...
    R0, B0, aMinor, Ip)
%EVALUATE_GLOBAL_PARAMETERS_LOCAL
% Compute volume, beta values and internal inductance.

    mu0 = 4*pi*1.e-7;

    psiNTable = psiNTable(:);
    pressureTable = pressureTable(:);
    psi = psi(:);

    if numel(psiNTable) ~= numel(pressureTable)
        error('GS:Finalize:PressureProfileSize', ...
            'psiN and pressure profiles must have equal lengths.');
    end

    if any(~isfinite(pressureTable))
        error('GS:Finalize:InvalidPressure', ...
            'Pressure contains NaN or Inf.');
    end

    pressureScale = ...
        max([max(abs(pressureTable)), 1]);

    if min(pressureTable) < -1.e-10*pressureScale
        error('GS:Finalize:NegativePressure', ...
            ['Pressure is negative. Check whether pedge is ', ...
             'defined as the boundary pressure.']);
    end

    if ~isfinite(R0) || R0 <= 0
        error('GS:Finalize:InvalidR0', ...
            'R0 must be positive.');
    end

    if ~isfinite(B0) || B0 == 0
        error('GS:Finalize:InvalidB0', ...
            'B0 must be finite and nonzero.');
    end

    currentScale = max(1, abs(Ip));

    if ~isfinite(Ip) || ...
            abs(Ip) <= eps(currentScale)
        error('GS:Finalize:InvalidCurrent', ...
            'Ip is zero or invalid.');
    end

    L1 = 1-quad.xi-quad.eta;
    L2 = quad.xi;
    L3 = quad.eta;

    plasmaVolume = 0;
    pressureIntegral = 0;
    Bp2VolumeIntegral = 0;

    nElement = size(P3.elements, 1);

    for e = 1:nElement
        vertexIDs = P3.P1elements(e, :);
        Ve = P3.P1points(vertexIDs, :);

        R1 = Ve(1, 1);
        Z1 = Ve(1, 2);

        R2 = Ve(2, 1);
        Z2 = Ve(2, 2);

        R3 = Ve(3, 1);
        Z3 = Ve(3, 2);

        J = [ ...
            R2-R1, R3-R1;
            Z2-Z1, Z3-Z1];

        detJ = det(J);

        if detJ <= 0
            error('GS:Finalize:InvalidElement', ...
                'Invalid element %d.', e);
        end

        Rq = L1*R1 + L2*R2 + L3*R3;

        if any(Rq <= 0)
            error('GS:Finalize:InvalidRadius', ...
                'Encountered R <= 0.');
        end

        solutionIDs = P3.elements(e, :);
        psiElement = psi(solutionIDs);

        % Flux at quadrature points
        psiQ = quad.N*psiElement;
        psiNQ = (psiQ-psiAxis)/dpsi;

        if any(psiNQ < -1.e-3) || ...
                any(psiNQ > 1+1.e-3)
            error('GS:Finalize:PsiNRange', ...
                ['Quadrature-point psiN extends too far ', ...
                 'outside [0,1].']);
        end

        psiNForProfile = ...
            min(max(psiNQ, 0), 1);

        pressureQ = interp1( ...
            psiNTable, pressureTable, ...
            psiNForProfile, 'linear');

        % Physical P3 derivatives
        dNdR = ( ...
              J(2, 2)*quad.dNdxi ...
            - J(2, 1)*quad.dNdeta)/detJ;

        dNdZ = ( ...
            - J(1, 2)*quad.dNdxi ...
            + J(1, 1)*quad.dNdeta)/detJ;

        dpsi_dR = dNdR*psiElement;
        dpsi_dZ = dNdZ*psiElement;

        BpSquared = ...
            (dpsi_dR.^2 + dpsi_dZ.^2)./Rq.^2;

        dA = quad.w*detJ;
        dV = 2*pi*Rq.*dA;

        plasmaVolume = ...
            plasmaVolume + sum(dV);

        pressureIntegral = ...
            pressureIntegral + sum(pressureQ.*dV);

        Bp2VolumeIntegral = ...
            Bp2VolumeIntegral + sum(BpSquared.*dV);
    end

    pressureVolumeAverage = ...
        pressureIntegral/plasmaVolume;

    betaTor = ...
        2*mu0*pressureVolumeAverage/B0^2;

    betaPol = ...
        4*pressureIntegral/(R0*mu0*Ip^2);

    li3 = ...
        2*Bp2VolumeIntegral/(R0*mu0^2*Ip^2);

    poloidalMagneticEnergy = ...
        Bp2VolumeIntegral/(2*mu0);

    result = struct();

    result.volume = plasmaVolume;

    result.pressureIntegral = ...
        pressureIntegral;

    result.pressureVolumeAverage = ...
        pressureVolumeAverage;

    result.Bp2VolumeIntegral = ...
        Bp2VolumeIntegral;

    result.poloidalMagneticEnergy = ...
        poloidalMagneticEnergy;

    result.betaTor = betaTor;
    result.betaPol = betaPol;
    result.li3 = li3;

    result.betaN = ...
        (100*betaTor)*aMinor*abs(B0)/(abs(Ip)/1.e6);
end


function assert_required_fields( ...
    value, fieldNames, variableName)
%ASSERT_REQUIRED_FIELDS Check required structure fields.

    if ~isstruct(value) || ~isscalar(value)
        error('GS:Finalize:InvalidStructure', ...
            '%s must be a scalar structure.', variableName);
    end

    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};

        if ~isfield(value, fieldName)
            error('GS:Finalize:MissingField', ...
                '%s.%s is required.', ...
                variableName, fieldName);
        end
    end
end
