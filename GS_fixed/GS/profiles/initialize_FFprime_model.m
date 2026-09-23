function model = initialize_FFprime_model( ...
    problem, profileInput, pprime, qTarget, options)
%INITIALIZE_FFPRIME_MODEL Initialize the mode-3 direct FF' predictor.
%
%   model = initialize_FFprime_model( ...
%       problem,profileInput,pprime,qTarget,options)
%
% The default initialization uses a one-dimensional cylindrical inverse
% estimate to determine the amplitude of the leading basis 1-psiN.  This
% estimate is only a preconditioner: all accepted updates in the mode-3
% solver are evaluated with the full two-dimensional GS solve.

    narginchk(4, 5);

    if nargin < 5 || isempty(options)
        options = struct();
    end

    assert_required_fields(problem, {'profile','dim'}, 'problem');
    assert_required_fields(problem.profile, {'psiN'}, 'problem.profile');
    assert_required_fields(problem.dim, {'R0','B0','aMinor'}, 'problem.dim');

    if ~isstruct(profileInput) || ~isscalar(profileInput)
        error('GS:FFprimeInit:ProfileInput', ...
            'profileInput must be a scalar structure.');
    end

    psiN = problem.profile.psiN(:);
    pprime = pprime(:);

    if numel(pprime) ~= numel(psiN) || any(~isfinite(pprime))
        error('GS:FFprimeInit:Pprime', ...
            'pprime must be finite and sampled on problem.profile.psiN.');
    end

    nBasis = get_option(options, 'nBasis', 3);
    enforceZeroEdge = logical(get_option( ...
        options, 'enforceZeroEdge', true));

    validateattributes(nBasis, {'numeric'}, ...
        {'real','finite','scalar','integer','positive','<=',12});

    [Fboundary, signF, normalization] = ...
        parse_fixed_boundary_normalization(problem, profileInput);

    basis = construct_FFprime_basis( ...
        psiN, nBasis, enforceZeroEdge);

    initial = get_initial_structure(profileInput);
    initializationType = get_initialization_type(initial);
    coefficients = zeros(nBasis,1);
    sourceFFprime = zeros(size(psiN));
    cylinderInfo = struct();

    switch initializationType
        case 'constant'
            % A constant G profile is exactly FF'=0.

        case 'ffprime'
            if ~isfield(initial,'FFprime') || isempty(initial.FFprime)
                error('GS:FFprimeInit:MissingFFprime', ...
                    'profileInput.initial.FFprime is required.');
            end

            initialFFprime = initial.FFprime(:);
            initialPsiN = get_initial_psiN( ...
                initial, psiN, numel(initialFFprime));

            sourceFFprime = interp1( ...
                initialPsiN, initialFFprime, psiN, 'pchip');

            coefficients = regularized_projection( ...
                basis, sourceFFprime, 1.0e-12);

        case 'cylindrical'
            qValues = evaluate_q_target(qTarget, psiN);

            [leadingCoefficient, sourceFFprime, cylinderInfo] = ...
                estimate_cylindrical_leading_coefficient( ...
                    psiN, pprime, qValues, Fboundary, ...
                    problem.dim.R0, problem.dim.aMinor, options);

            coefficients(1) = leadingCoefficient;

        otherwise
            error('GS:FFprimeInit:InitializationType', ...
                'Unknown initialization type: %s.', initializationType);
    end

    if any(~isfinite(coefficients))
        error('GS:FFprimeInit:InvalidCoefficients', ...
            'The initialized FFprime coefficients are nonfinite.');
    end

    model = struct();
    model.representation = 'direct-FFprime-reduced-basis';
    model.psiN = psiN;
    model.nBasis = nBasis;
    model.enforceZeroEdge = enforceZeroEdge;
    model.coefficients = coefficients;
    model.Fboundary = Fboundary;
    model.signF = signF;
    model.normalization = normalization;

    model.initialization = struct();
    model.initialization.type = initializationType;
    model.initialization.sourceFFprime = sourceFFprime;
    model.initialization.cylindrical = cylinderInfo;

    % Validate positivity of the reconstructed initial profile here.
    evaluate_FFprime_model(model, psiN);
end


function [coefficient, FFprimeCylinder, info] = ...
    estimate_cylindrical_leading_coefficient( ...
        psiN, pprime, qTarget, Fboundary, R0, aMinor, options)
%ESTIMATE_CYLINDRICAL_LEADING_COEFFICIENT
% Construct a cheap circular-cylinder estimate and project it onto 1-s.

    qTarget = qTarget(:);

    if numel(qTarget) ~= numel(psiN) || ...
            any(~isfinite(qTarget)) || any(qTarget == 0)
        error('GS:FFprimeInit:InvalidQTarget', ...
            'The cylindrical initializer requires finite nonzero q.');
    end

    nIteration = get_option(options, 'cylinderIterations', 6);
    relaxation = get_option(options, 'cylinderRelaxation', 0.5);
    maximumFraction = get_option(options, 'cylinderMaximumGDrop', 0.8);

    validateattributes(nIteration, {'numeric'}, ...
        {'real','finite','scalar','integer','positive'});
    validateattributes(relaxation, {'numeric'}, ...
        {'real','finite','scalar','>',0,'<=',1});
    validateattributes(maximumFraction, {'numeric'}, ...
        {'real','finite','scalar','>',0,'<',1});

    qMagnitude = abs(qTarget);
    FboundaryMagnitude = abs(Fboundary);
    F = FboundaryMagnitude*ones(size(psiN));
    coefficient = 0;
    mu0 = 4*pi*1.0e-7;
    leadingBasis = 1-psiN;
    fitMask = psiN <= 0.9;

    coefficientHistory = nan(nIteration,1);
    fluxSpanHistory = nan(nIteration,1);
    FFprimeCylinder = zeros(size(psiN));

    for iteration = 1:nIteration
        qOverFIntegral = trapz(psiN, qMagnitude./F);

        if ~isfinite(qOverFIntegral) || qOverFIntegral <= 0
            error('GS:FFprimeInit:CylinderIntegral', ...
                'Failed to estimate the cylindrical flux span.');
        end

        fluxSpanMagnitude = ...
            aMinor^2/(2*R0*qOverFIntegral);

        radiusSquared = 2*R0*fluxSpanMagnitude ...
            *cumtrapz(psiN, qMagnitude./F);

        if radiusSquared(end) <= 0 || ...
                any(~isfinite(radiusSquared))
            error('GS:FFprimeInit:CylinderRadius', ...
                'Failed to reconstruct the cylindrical radius.');
        end

        radiusSquared = ...
            radiusSquared*(aMinor^2/radiusSquared(end));

        cylindricalFluxFunction = ...
            radiusSquared.*F./(R0^2*qMagnitude);

        derivativePP = differentiate_scalar_pp( ...
            pchip(psiN,cylindricalFluxFunction));

        derivativeFluxFunction = ...
            ppval(derivativePP,psiN).';

        % The fixed-boundary solver normally uses a maximum-type magnetic
        % axis, so dpsi=psi_b-psi_a<0.  Ampere's law and the cylindrical
        % safety-factor relation then give the following normalized source.
        FFprimeCylinder = ...
            -F./qMagnitude.*derivativeFluxFunction ...
            -mu0*R0^2*pprime;

        numerator = leadingBasis(fitMask).'* ...
            FFprimeCylinder(fitMask);
        denominator = leadingBasis(fitMask).'* ...
            leadingBasis(fitMask);

        candidate = numerator/max(denominator,eps);

        maximumPositive = ...
            maximumFraction*FboundaryMagnitude^2;
        maximumMagnitude = FboundaryMagnitude^2;

        candidate = min(candidate, maximumPositive);
        candidate = max(candidate, -maximumMagnitude);

        coefficient = ...
            (1-relaxation)*coefficient + relaxation*candidate;

        trialFFprime = coefficient*leadingBasis;
        reconstructed = reconstruct_G_from_FFprime( ...
            psiN, trialFFprime, Fboundary);

        F = abs(reconstructed.F);
        coefficientHistory(iteration) = coefficient;
        fluxSpanHistory(iteration) = fluxSpanMagnitude;
    end

    info = struct();
    info.leadingCoefficient = coefficient;
    info.coefficientHistory = coefficientHistory;
    info.fluxSpanMagnitudeHistory = fluxSpanHistory;
    info.finalCylinderFFprime = FFprimeCylinder;
end


function coefficients = regularized_projection(basis, values, ridge)
%REGULARIZED_PROJECTION Project samples onto the reduced basis.

    normalMatrix = basis.'*basis + ridge*eye(size(basis,2));
    coefficients = normalMatrix\(basis.'*values(:));
end


function [Fboundary, signF, normalization] = ...
    parse_fixed_boundary_normalization(problem, profileInput)
%PARSE_FIXED_BOUNDARY_NORMALIZATION Parse the mode-3 F normalization.

    defaultFboundary = problem.dim.R0*problem.dim.B0;

    if isfield(profileInput,'normalization') && ...
            ~isempty(profileInput.normalization)
        normalizationInput = profileInput.normalization;
    else
        normalizationInput = struct();
    end

    if ~isstruct(normalizationInput) || ~isscalar(normalizationInput)
        error('GS:FFprimeInit:Normalization', ...
            'profileInput.normalization must be a scalar structure.');
    end

    mode = get_option(normalizationInput, 'mode', 'fixed_Fb');

    if isstring(mode)
        mode = char(mode);
    end

    mode = lower(strrep(strtrim(mode),'-','_'));

    if strcmp(mode,'fixedfb')
        mode = 'fixed_fb';
    end

    if ~strcmp(mode,'fixed_fb')
        error('GS:FFprimeInit:NormalizationMode', ...
            ['The direct-FFprime mode-3 solver currently requires ', ...
             'normalization.mode=''fixed_Fb''.']);
    end

    Fboundary = get_option( ...
        normalizationInput, 'Fb', defaultFboundary);

    validateattributes(Fboundary, {'numeric'}, ...
        {'real','finite','scalar','nonzero'});

    if isfield(profileInput,'FSign') && ...
            ~isempty(profileInput.FSign)
        signF = profileInput.FSign;
    else
        signF = sign(Fboundary);
    end

    if ~(isscalar(signF) && isfinite(signF) && ...
            (signF == 1 || signF == -1))
        error('GS:FFprimeInit:FSign', ...
            'profileInput.FSign must be +1 or -1.');
    end

    if sign(Fboundary) ~= signF
        error('GS:FFprimeInit:FSignMismatch', ...
            'normalization.Fb and profileInput.FSign are inconsistent.');
    end

    normalization = struct( ...
        'mode','fixed_fb', ...
        'Fb',Fboundary, ...
        'signF',signF);
end


function initial = get_initial_structure(profileInput)
%GET_INITIAL_STRUCTURE Return optional initialization data.

    if isfield(profileInput,'initial') && ...
            ~isempty(profileInput.initial)
        initial = profileInput.initial;
    else
        initial = struct();
    end

    if ~isstruct(initial) || ~isscalar(initial)
        error('GS:FFprimeInit:InitialStructure', ...
            'profileInput.initial must be a scalar structure.');
    end
end


function initializationType = get_initialization_type(initial)
%GET_INITIALIZATION_TYPE Resolve the initialization strategy.

    initializationType = get_option(initial, 'type', 'cylindrical');

    if isstring(initializationType)
        initializationType = char(initializationType);
    end

    if ~ischar(initializationType) || ~isrow(initializationType)
        error('GS:FFprimeInit:InitializationType', ...
            'profileInput.initial.type must be scalar text.');
    end

    initializationType = lower(strtrim(initializationType));

    if strcmp(initializationType,'zero') || ...
            strcmp(initializationType,'constant_g')
        initializationType = 'constant';
    end

    if ~ismember(initializationType, ...
            {'constant','ffprime','cylindrical'})
        error('GS:FFprimeInit:InitializationType', ...
            ['profileInput.initial.type must be constant, ', ...
             'FFprime, or cylindrical.']);
    end
end


function psiN = get_initial_psiN(initial, commonPsiN, dataLength)
%GET_INITIAL_PSIN Resolve an optional input profile grid.

    if isfield(initial,'psiN') && ~isempty(initial.psiN)
        psiN = initial.psiN(:);
    else
        if dataLength ~= numel(commonPsiN)
            error('GS:FFprimeInit:MissingInitialPsiN', ...
                'initial.psiN is required for a nonstandard profile size.');
        end
        psiN = commonPsiN;
    end

    if numel(psiN) ~= dataLength || any(~isfinite(psiN)) || ...
            any(diff(psiN) <= 0)
        error('GS:FFprimeInit:InitialProfileGrid', ...
            'The initial profile grid is invalid.');
    end
end


function values = evaluate_q_target(qTarget, psiN)
%EVALUATE_Q_TARGET Evaluate supported prescribed-q representations.

    if isa(qTarget,'function_handle')
        values = qTarget(psiN);
    elseif isstruct(qTarget)
        if ~isfield(qTarget,'psiN') || ~isfield(qTarget,'q')
            error('GS:FFprimeInit:QTargetStructure', ...
                'qTarget must contain psiN and q.');
        end
        values = interp1( ...
            qTarget.psiN(:), qTarget.q(:), psiN, 'pchip');
    elseif isscalar(qTarget)
        values = qTarget*ones(size(psiN));
    elseif isnumeric(qTarget) && numel(qTarget) == numel(psiN)
        values = qTarget(:);
    else
        error('GS:FFprimeInit:QTarget', ...
            'Unsupported qTarget representation.');
    end

    values = values(:);

    if numel(values) ~= numel(psiN) || ...
            any(~isfinite(values)) || any(values == 0)
        error('GS:FFprimeInit:QTarget', ...
            'The evaluated q target must be finite and nonzero.');
    end
end


function derivativePP = differentiate_scalar_pp(pp)
%DIFFERENTIATE_SCALAR_PP Differentiate a scalar pp representation.

    [breaks, coefficients, nPiece, order, dimension] = unmkpp(pp);

    if dimension ~= 1
        error('GS:FFprimeInit:PPDimension', ...
            'Only scalar piecewise polynomials are supported.');
    end

    if order <= 1
        derivativePP = mkpp(breaks,zeros(nPiece,1));
        return
    end

    powers = order-1:-1:1;
    derivativePP = mkpp( ...
        breaks, coefficients(:,1:order-1).*powers);
end


function value = get_option(options, fieldName, defaultValue)
%GET_OPTION Read an optional scalar-structure field.

    if isfield(options,fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
    else
        value = defaultValue;
    end
end


function assert_required_fields(value, requiredFields, valueName)
%ASSERT_REQUIRED_FIELDS Validate a scalar structure interface.

    if ~isstruct(value) || ~isscalar(value)
        error('GS:FFprimeInit:InvalidStructure', ...
            '%s must be a scalar structure.', valueName);
    end

    for k = 1:numel(requiredFields)
        if ~isfield(value, requiredFields{k})
            error('GS:FFprimeInit:MissingField', ...
                '%s.%s is required.', valueName, requiredFields{k});
        end
    end
end
