function point = evaluate_salpha_point(state, shat, alpha, opts)
%EVALUATE_SALPHA_POINT Solve one (shat,alpha) point over ballooning phase.
%
%   point = EVALUATE_SALPHA_POINT(state, shat, alpha)
%
% state is returned by prepare_salpha_equilibrium or
% prepare_fem_salpha_equilibrium.  point.lambdaMax is the largest
% ballooning eigenvalue over the requested Theta0 phases.

    arguments
        state (1,1) struct
        shat (1,1) double
        alpha (1,1) double
        opts.Theta0 double = []
        opts.ThetaBoundary (1,1) double = 5*pi
        opts.NBalloonPoints (1,1) double = 151
        opts.EnforcePeriodicity (1,1) logical = true
    end

    validateattributes(shat, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'shat');
    validateattributes(alpha, {'numeric'}, ...
        {'real', 'finite', 'scalar'}, mfilename, 'alpha');
    validateattributes(opts.ThetaBoundary, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'ThetaBoundary');
    validateattributes(opts.NBalloonPoints, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 21}, ...
        mfilename, 'NBalloonPoints');

    if mod(opts.NBalloonPoints, 2) == 0
        error('salpha:Point:EvenBalloonGrid', ...
            'NBalloonPoints must be odd so theta=0 is represented.');
    end

    if isempty(opts.Theta0)
        nTheta0 = 12;
        theta0 = 2*pi*((0:nTheta0-1)-floor(nTheta0/2))/nTheta0;
    else
        theta0 = opts.Theta0(:).';
        if any(~isfinite(theta0)) || ~isreal(theta0) || ...
                numel(unique(theta0)) ~= numel(theta0)
            error('salpha:Point:InvalidTheta0', ...
                'Theta0 must contain distinct, finite real values.');
        end
    end

    requiredStateFields = { ...
        'localEquilibrium', 'integrals', 'surface', ...
        'FFprime', 'normalization'};
    assert_fields(state, requiredStateFields);

    volume = state.integrals.V;
    Vprime = state.integrals.Vprime;
    q = state.surface.q;
    Raxis = state.surface.Raxis;

    if ~all(isfinite([volume, Vprime, q, Raxis])) || ...
            volume <= 0 || Vprime <= 0 || q == 0 || Raxis <= 0
        error('salpha:Point:InvalidSurfaceIntegral', ...
            'V, Vprime, q, and Raxis must define a valid flux surface.');
    end

    mu0 = 4*pi*1.e-7;
    qprime = shat*(q*Vprime/(2*volume));
    pprime = -alpha*(4*pi^2/(2*Vprime)) ...
        *sqrt(2*pi^2*Raxis/volume)/mu0;

    [~, ~, thetaDiagnostic] = evaluate_theta_psi( ...
        state.localEquilibrium, state.FFprime, pprime, qprime, ...
        EnforcePeriodicity=opts.EnforcePeriodicity);

    thetaGrid = linspace( ...
        -opts.ThetaBoundary, opts.ThetaBoundary, ...
        opts.NBalloonPoints).';
    lambdaTheta0 = nan(size(theta0));

    for thetaIndex = 1:numel(theta0)
        coefficient = evaluate_gcf( ...
            thetaDiagnostic, theta0(thetaIndex), ...
            theta_bnd=opts.ThetaBoundary, ...
            a_N=state.normalization.aN, ...
            B_N=state.normalization.BN);

        g = interp1( ...
            coefficient.theta, coefficient.g, thetaGrid, 'pchip');
        c = interp1( ...
            coefficient.theta, coefficient.c, thetaGrid, 'pchip');
        f = interp1( ...
            coefficient.theta, coefficient.f, thetaGrid, 'pchip');

        [Gmatrix, Cmatrix, Fmatrix, thetaDof] = ...
            construct_matrix(thetaGrid, g, c, f);
        lambda = calculate_eigenvalue( ...
            thetaDof, Gmatrix, Cmatrix, Fmatrix, 'dirichlet');
        lambdaTheta0(thetaIndex) = lambda(1);
    end

    [lambdaMaximum, maximumIndex] = max(lambdaTheta0);

    point = struct();
    point.shat = shat;
    point.alpha = alpha;
    point.theta0 = theta0;
    point.lambdaTheta0 = lambdaTheta0;
    point.lambdaMax = lambdaMaximum;
    point.bestTheta0 = theta0(maximumIndex);
    point.qprime = qprime;
    point.pprime = pprime;
    point.FFprimeUsed = thetaDiagnostic.FFprime_used;
    point.periodicityClosure = thetaDiagnostic.closure_used;
end


function assert_fields(value, names)
    if ~isstruct(value) || ~isscalar(value)
        error('salpha:Point:InvalidState', ...
            'state must be a scalar structure.');
    end

    for index = 1:numel(names)
        if ~isfield(value, names{index})
            error('salpha:Point:MissingField', ...
                'state.%s is required.', names{index});
        end
    end
end
