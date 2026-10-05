function result = scan_salpha_diagram(state, opts)
%SCAN_SALPHA_DIAGRAM Scan ballooning stability at fixed surface geometry.
%
%   result = SCAN_SALPHA_DIAGRAM(state)
%
% state is produced by prepare_salpha_equilibrium or
% prepare_fem_salpha_equilibrium.  The surface shape, Bp, F, and q are
% held fixed while qprime and pprime are chosen from the requested
% magnetic-shear and pressure-gradient coordinates,
%
%   shat = 2*V*qprime/(q*Vprime),
%   alpha = -mu0*pprime*(2*Vprime)/(4*pi^2)
%           *sqrt(V/(2*pi^2*Raxis)).
%
% The existing evaluate_theta_psi, evaluate_gcf, construct_matrix, and
% calculate_eigenvalue routines are used without modification.  At each
% (shat,alpha), result.lambda_max is the largest eigenvalue over Theta0.
% Positive lambda_max denotes ballooning instability.

    arguments
        state (1,1) struct
        opts.Shat double = linspace(0, 5, 31)
        opts.Alpha double = linspace(0, 6, 31)
        opts.Theta0 double = []
        opts.ThetaBoundary (1,1) double = 5*pi
        opts.NBalloonPoints (1,1) double = 151
        opts.EnforcePeriodicity (1,1) logical = true
        opts.Verbose (1,1) logical = true
    end

    requiredStateFields = { ...
        'localEquilibrium', 'integrals', 'surface', ...
        'FFprime', 'normalization', 'targetPsiN'};
    assert_fields(state, requiredStateFields, 'state');

    requiredIntegralFields = {'V', 'Vprime'};
    assert_fields(state.integrals, requiredIntegralFields, ...
        'state.integrals');

    requiredSurfaceFields = {'q', 'Raxis'};
    assert_fields(state.surface, requiredSurfaceFields, ...
        'state.surface');

    requiredNormalizationFields = {'aN', 'BN'};
    assert_fields(state.normalization, requiredNormalizationFields, ...
        'state.normalization');

    shat = validate_scan_axis(opts.Shat, 'Shat');
    alpha = validate_scan_axis(opts.Alpha, 'Alpha');

    if isempty(opts.Theta0)
        nTheta0Default = 12;
        theta0 = 2*pi*((0:nTheta0Default-1) ...
            - floor(nTheta0Default/2))/nTheta0Default;
    else
        theta0 = opts.Theta0(:).';
        if any(~isfinite(theta0)) || ~isreal(theta0) || ...
                numel(unique(theta0)) ~= numel(theta0)
            error('salpha:Scan:InvalidTheta0', ...
                'Theta0 must contain distinct, finite real values.');
        end
    end

    validateattributes(opts.ThetaBoundary, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'ThetaBoundary');

    validateattributes(opts.NBalloonPoints, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 21}, ...
        mfilename, 'NBalloonPoints');

    if mod(opts.NBalloonPoints, 2) == 0
        error('salpha:Scan:EvenBalloonGrid', ...
            'NBalloonPoints must be odd so theta=0 is represented.');
    end

    volume = state.integrals.V;
    Vprime = state.integrals.Vprime;
    q = state.surface.q;
    Raxis = state.surface.Raxis;

    if ~all(isfinite([volume, Vprime, q, Raxis])) || ...
            volume <= 0 || Vprime <= 0 || q == 0 || Raxis <= 0
        error('salpha:Scan:InvalidSurfaceIntegral', ...
            'V, Vprime, q, and Raxis must define a valid flux surface.');
    end

    mu0 = 4*pi*1.e-7;
    qprimeScan = shat.*(q*Vprime/(2*volume));
    pprimeScan = -alpha.*(4*pi^2/(2*Vprime)) ...
        .*sqrt(2*pi^2*Raxis/volume)/mu0;

    thetaGrid = linspace( ...
        -opts.ThetaBoundary, opts.ThetaBoundary, ...
        opts.NBalloonPoints).';

    nShat = numel(shat);
    nAlpha = numel(alpha);
    nTheta0 = numel(theta0);

    lambdaMaximum = nan(nShat, nAlpha);
    bestTheta0 = nan(nShat, nAlpha);
    lambdaTheta0 = nan(nShat, nAlpha, nTheta0);
    failMessage = strings(nShat, nAlpha, nTheta0);
    FFprimeUsed = nan(nShat, nAlpha);
    periodicityClosure = nan(nShat, nAlpha);

    if opts.Verbose
        fprintf(['Starting s-alpha scan: %d shat x %d alpha x ', ...
            '%d theta0 = %d solves.\n'], ...
            nShat, nAlpha, nTheta0, nShat*nAlpha*nTheta0);
    end

    scanTimer = tic;
    for shatIndex = 1:nShat
        rowTimer = tic;
        qprime = qprimeScan(shatIndex);

        for alphaIndex = 1:nAlpha
            pprime = pprimeScan(alphaIndex);

            try
                [~, ~, thetaDiagnostic] = evaluate_theta_psi( ...
                    state.localEquilibrium, ...
                    state.FFprime, pprime, qprime, ...
                    EnforcePeriodicity=opts.EnforcePeriodicity);

                FFprimeUsed(shatIndex, alphaIndex) = ...
                    thetaDiagnostic.FFprime_used;
                periodicityClosure(shatIndex, alphaIndex) = ...
                    thetaDiagnostic.closure_used;
            catch exception
                for thetaIndex = 1:nTheta0
                    failMessage( ...
                        shatIndex, alphaIndex, thetaIndex) = ...
                        string(exception.message);
                end
                continue;
            end

            for thetaIndex = 1:nTheta0
                try
                    coefficient = evaluate_gcf( ...
                        thetaDiagnostic, theta0(thetaIndex), ...
                        theta_bnd=opts.ThetaBoundary, ...
                        a_N=state.normalization.aN, ...
                        B_N=state.normalization.BN);

                    g = interp1( ...
                        coefficient.theta, coefficient.g, ...
                        thetaGrid, 'pchip');
                    c = interp1( ...
                        coefficient.theta, coefficient.c, ...
                        thetaGrid, 'pchip');
                    f = interp1( ...
                        coefficient.theta, coefficient.f, ...
                        thetaGrid, 'pchip');

                    [Gmatrix, Cmatrix, Fmatrix, thetaDof] = ...
                        construct_matrix(thetaGrid, g, c, f);

                    lambda = calculate_eigenvalue( ...
                        thetaDof, Gmatrix, Cmatrix, Fmatrix, ...
                        'dirichlet');

                    lambdaTheta0( ...
                        shatIndex, alphaIndex, thetaIndex) = ...
                        lambda(1);
                catch exception
                    failMessage( ...
                        shatIndex, alphaIndex, thetaIndex) = ...
                        string(exception.message);
                end
            end

            values = reshape( ...
                lambdaTheta0(shatIndex, alphaIndex, :), 1, []);
            [lambdaMaximum(shatIndex, alphaIndex), maximumIndex] = ...
                max(values, [], 'omitnan');

            if isfinite(lambdaMaximum(shatIndex, alphaIndex))
                bestTheta0(shatIndex, alphaIndex) = ...
                    theta0(maximumIndex);
            end
        end

        if opts.Verbose
            elapsed = toc(scanTimer);
            eta = elapsed/shatIndex*(nShat-shatIndex);
            fprintf(['  shat row %d/%d: shat=%+.5g, ', ...
                'row %.1f s, ETA %.1f min\n'], ...
                shatIndex, nShat, shat(shatIndex), ...
                toc(rowTimer), eta/60);
        end
    end

    elapsedSeconds = toc(scanTimer);
    numberFailed = nnz(strlength(failMessage) > 0);
    finiteLambda = lambdaMaximum(isfinite(lambdaMaximum));

    if isempty(finiteLambda)
        lambdaRange = [NaN, NaN];
        hasMarginalCrossing = false;
    else
        lambdaRange = [min(finiteLambda), max(finiteLambda)];
        hasMarginalCrossing = ...
            lambdaRange(1) <= 0 && lambdaRange(2) >= 0;
    end

    result = struct();
    if isfield(state.localEquilibrium, 'model')
        result.model = state.localEquilibrium.model;
    elseif isfield(state, 'model')
        result.model = state.model;
    else
        result.model = 'unspecified-local-equilibrium';
    end
    result.targetPsiN = state.targetPsiN;

    result.shat = shat;
    result.alpha = alpha;
    result.theta0 = theta0;
    result.thetaGrid = thetaGrid;

    result.lambdaMax = lambdaMaximum;
    result.lambdaTheta0 = lambdaTheta0;
    result.bestTheta0 = bestTheta0;
    result.failMessage = failMessage;

    result.qprime = qprimeScan;
    result.pprime = pprimeScan;
    result.FFprimeInput = state.FFprime;
    result.FFprimeUsed = FFprimeUsed;
    result.periodicityClosure = periodicityClosure;

    result.equilibriumShat = state.shat;
    result.equilibriumAlpha = state.alpha;
    result.volume = volume;
    result.Vprime = Vprime;
    result.q = q;

    result.lambdaRange = lambdaRange;
    result.hasMarginalCrossing = hasMarginalCrossing;
    result.numberFailed = numberFailed;
    result.scanElapsedSeconds = elapsedSeconds;
    result.settings = struct( ...
        'thetaBoundary', opts.ThetaBoundary, ...
        'nBalloonPoints', opts.NBalloonPoints, ...
        'periodicityEnforced', opts.EnforcePeriodicity);

    if opts.Verbose
        fprintf(['s-alpha scan finished in %.1f s: ', ...
            '%d failed solves, lambda range [%+.5e,%+.5e].\n'], ...
            elapsedSeconds, numberFailed, ...
            lambdaRange(1), lambdaRange(2));
    end
end


function axisValue = validate_scan_axis(value, name)
    if ~isvector(value) || numel(value) < 2 || ...
            ~isreal(value) || any(~isfinite(value))
        error('salpha:Scan:InvalidAxis', ...
            '%s must contain at least two finite real values.', name);
    end

    axisValue = value(:).';
    if any(diff(axisValue) <= 0)
        error('salpha:Scan:UnorderedAxis', ...
            '%s must be strictly increasing.', name);
    end
end


function assert_fields(value, names, variableName)
    if ~isstruct(value) || ~isscalar(value)
        error('salpha:Scan:InvalidStructure', ...
            '%s must be a scalar structure.', variableName);
    end

    for index = 1:numel(names)
        name = names{index};
        if ~isfield(value, name)
            error('salpha:Scan:MissingField', ...
                '%s.%s is required.', variableName, name);
        end
    end
end
