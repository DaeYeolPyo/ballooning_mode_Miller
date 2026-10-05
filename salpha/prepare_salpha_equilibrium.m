function state = prepare_salpha_equilibrium( ...
    eq, eqfunc, targetPsiN, opts)
%PREPARE_SALPHA_EQUILIBRIUM Build one local-equilibrium ballooning state.
%
%   state = PREPARE_SALPHA_EQUILIBRIUM(eq, eqfunc, targetPsiN)
%
% This function collects the equilibrium-dependent setup that precedes an
% s-alpha scan.  It accepts the same eq/eqfunc interface used by the GEQDSK
% path, so callers may supply either read_geqdsk/build_interpolants output
% or the FEM adapter output from build_ballooning_adapter.
%
% Model="Miller" uses the existing analytic local fit.  Model="Direct"
% constructs the Mercier-Luc surface quantities directly from the extracted
% contour and eqfunc gradient, without assuming a Miller shape.  The
% returned thetaDiag can be passed unchanged to evaluate_gcf,
% construct_matrix, and calculate_eigenvalue.

    arguments
        eq (1,1) struct
        eqfunc (1,1) struct
        targetPsiN (1,1) double
        opts.Model (1,1) string = "Miller"
        opts.NSurfacePoints (1,1) double = 1024
        opts.RadialFitHalfWidth (1,1) double = 3.e-2
        opts.NRadialFit (1,1) double = 9
        opts.EnforcePeriodicity (1,1) logical = true
        opts.FourierModes (1,1) double = 32
        opts.ProjectionIterations (1,1) double = 2
        opts.EnforceDirectQ (1,1) logical = true
    end

    validateattributes(targetPsiN, {'numeric'}, ...
        {'real', 'finite', 'scalar', '>', 0, '<', 1}, ...
        mfilename, 'targetPsiN');

    validateattributes(opts.NSurfacePoints, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 32}, ...
        mfilename, 'NSurfacePoints');

    validateattributes(opts.RadialFitHalfWidth, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, ...
        mfilename, 'RadialFitHalfWidth');

    validateattributes(opts.NRadialFit, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 3}, ...
        mfilename, 'NRadialFit');

    if mod(opts.NRadialFit, 2) == 0
        error('salpha:PrepareEquilibrium:EvenRadialFit', ...
            'NRadialFit must be odd.');
    end

    model = lower(strtrim(opts.Model));
    if ~ismember(model, ["miller", "direct", "flux-surface"])
        error('salpha:PrepareEquilibrium:UnsupportedModel', ...
            ['Supported models are "Miller" and "Direct" ', ...
             '(alias "flux-surface"). Requested model: %s.'], ...
            opts.Model);
    end

    surface = extract_flux_surface( ...
        eq, eqfunc, targetPsiN, ...
        Npoints=int32(opts.NSurfacePoints));

    integrals = compute_contour_integrals( ...
        eq, eqfunc, surface, isNormalized=false);

    switch model
        case "miller"
            [parameters, boundary] = fit_Miller( ...
                eq, eqfunc, targetPsiN, ...
                NTheta=int32(opts.NSurfacePoints), ...
                dpsi=opts.RadialFitHalfWidth, ...
                NRadialFit=int32(opts.NRadialFit));

            localEquilibrium = Miller_Mercier_Luc( ...
                surface, parameters, boundary);

        case {"direct", "flux-surface"}
            localEquilibrium = FluxSurface_Mercier_Luc( ...
                surface, eqfunc, ...
                FourierModes=opts.FourierModes, ...
                ProjectionIterations=opts.ProjectionIterations, ...
                EnforceQ=opts.EnforceDirectQ);

            parameters = struct();
            parameters.FourierModes = opts.FourierModes;
            parameters.ProjectionIterations = ...
                opts.ProjectionIterations;

            boundary = struct();
            boundary.R = localEquilibrium.R;
            boundary.Z = localEquilibrium.Z;
            boundary.Rt = localEquilibrium.Rt;
            boundary.Zt = localEquilibrium.Zt;
            boundary.Rtt = localEquilibrium.Rtt;
            boundary.Ztt = localEquilibrium.Ztt;
            boundary.theta = localEquilibrium.theta_geo;
    end

    FFprime = eqfunc.FFprime(targetPsiN);
    pprime = eqfunc.pprime(targetPsiN);
    qprime = eqfunc.qprime(targetPsiN);

    [Gamma, thetaPsi, thetaDiag] = evaluate_theta_psi( ...
        localEquilibrium, FFprime, pprime, qprime, ...
        EnforcePeriodicity=opts.EnforcePeriodicity);

    state = struct();
    state.model = char(model);
    state.targetPsiN = targetPsiN;

    state.surface = surface;
    state.integrals = integrals;

    state.parameters = parameters;
    state.boundary = boundary;
    state.localEquilibrium = localEquilibrium;

    state.FFprime = FFprime;
    state.pprime = pprime;
    state.qprime = qprime;

    state.shat = integrals.s_hat;
    state.alpha = integrals.alpha;

    state.Gamma = Gamma;
    state.thetaPsi = thetaPsi;
    state.thetaDiag = thetaDiag;

    state.normalization = struct();
    state.normalization.aN = eq.rmaxis;
    state.normalization.BN = abs(eq.bcentr);

    state.diagnostics = struct();
    state.diagnostics.qIntegralRelativeError = ...
        integrals.q_relerr;
    state.diagnostics.qLocalRelativeError = ...
        (localEquilibrium.q_check-surface.q) ...
        /max(abs(surface.q), eps);
    state.diagnostics.FFprimePeriodicityRelativeDifference = ...
        abs(thetaDiag.FFprime_periodic-FFprime) ...
        /max(abs(FFprime), 1);
    state.diagnostics.periodicityClosure = ...
        thetaDiag.closure_used;

    if isfield(localEquilibrium, 'diagnostics')
        state.diagnostics.directGeometry = ...
            localEquilibrium.diagnostics;
    end
end
