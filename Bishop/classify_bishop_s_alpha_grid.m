function grid = classify_bishop_s_alpha_grid(geom, varargin)
%CLASSIFY_BISHOP_S_ALPHA_GRID Classify stable/unstable on an (s, alpha) grid.
%
%   grid = CLASSIFY_BISHOP_S_ALPHA_GRID(geom)
%   grid = CLASSIFY_BISHOP_S_ALPHA_GRID(geom, 'Name', value, ...)
%
% For each supplied (s, alpha) point this function converts
%
%   s     -> Iprime
%   alpha -> pprimeEquation
%
% then evaluates the Bishop Chapter 3 ballooning operator and classifies
% the point by the sign of lambda_min.
%
% Definitions used by default
%   q ~= I*Cq,
%   Cq = (1/2pi) int dl/(R^2 Bp),
%
%   s = 2*(psi - psi_axis)/q * dq/dpsi,
%   The default implementation does not use dq/dpsi ~= Iprime*Cq directly
%   to set Iprime. Instead, it sets Iprime so that the one-period
%   ballooning secular increment Jperiod in Eq. (29) matches the requested
%   magnetic shear. This keeps s and alpha independent even though p'
%   also contributes to Jperiod.
%
%   alpha = AlphaFactor * pprimeEquation,
%   pprimeEquation = mu0 * dp/dpsi_GEQDSK,
%
%   AlphaFactor = -2*q^2*R0*<|grad psi|>/<B^2>.
%
% Options
%   'SValues'              : vector of s points. Default linspace(-2,2,9).
%   'AlphaValues'          : vector of alpha points. Default linspace(0,1,9).
%   'I0'                   : override I. Default geom.I0.
%   'Q0'                   : q used in labels/conversion. Default geom.q0
%                            if available, otherwise I0*Cq.
%   'ShearFactor'          : factor in s = factor*qprime/q. Default
%                            2*(psi_target - psi_axis).
%   'AlphaFactor'          : factor in alpha = factor*pprimeEquation.
%                            Default large-aspect-ratio estimate above.
%   'ShearMapping'         : 'jperiod' or 'qprime'. Default 'jperiod'.
%                            'jperiod' solves Iprime from the desired
%                            Eq. (29) Jperiod at every (s,alpha).
%                            'qprime' uses the simple Iprime=qprime/Cq.
%   'JShearFactor'         : Jperiod per unit s for ShearMapping='jperiod'.
%                            Default is calibrated from geom.shear0 and
%                            the baseline geom.Iprime0, geom.mu0_pprime0.
%   'NPeriodsEachSide'     : ballooning theta half-window. Default 5.
%   'ScanBallooningPhase'  : true scans theta_k/L0Index and keeps the most
%                            unstable eigenvalue, as required by Bishop's
%                            choice of l0 below Eq. (25). Default true.
%   'BallooningPhaseIndices': explicit L0Index values for theta_k scan.
%   'NumBallooningPhases'  : number of approximately uniform phases if
%                            ScanBallooningPhase=true and explicit indices
%                            are not supplied. Default 16.
%   'NumModes'             : modes requested from stability solve. Default 3.
%   'StabilityTolerance'   : lambda tolerance around zero. Default 1e-8.
%   'Plot'                 : true/false. Default true.
%   'Verbose'              : true/false. Default true.
%
% Output
%   grid.lambdaMin(i,j)    : min_l0 lowest eigenvalue at s_i, alpha_j
%   grid.phaseIndexMin(i,j): L0Index that produced lambdaMin
%   grid.phaseL0Min(i,j)   : l0 value that produced lambdaMin
%   grid.phaseThetaMin(i,j): 2*pi*l0/L for that l0
%   grid.status(i,j)       : "stable", "unstable", "marginal", or "failed"
%   grid.stable            : logical stable mask
%   grid.unstable          : logical unstable mask

    p = inputParser;
    addParameter(p, 'SValues', linspace(-2, 2, 9), @(x)isnumeric(x) && isvector(x));
    addParameter(p, 'AlphaValues', linspace(0, 1, 9), @(x)isnumeric(x) && isvector(x));
    addParameter(p, 'I0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Q0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'ShearFactor', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'AlphaFactor', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'ShearMapping', 'jperiod', @(x)ischar(x) || isstring(x));
    addParameter(p, 'JShearFactor', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'L0Index', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'ScanBallooningPhase', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'BallooningPhaseIndices', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'NumBallooningPhases', 16, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NPeriodsEachSide', 5, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NumModes', 3, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'StabilityTolerance', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Plot', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Verbose', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    constants = compute_constants(geom, opt);
    shearMapping = lower(string(opt.ShearMapping));
    if ~any(shearMapping == ["jperiod", "qprime"])
        error('classify_bishop_s_alpha_grid:BadShearMapping', ...
            'ShearMapping must be ''jperiod'' or ''qprime''.');
    end

    sValues = opt.SValues(:);
    alphaValues = opt.AlphaValues(:).';
    ns = numel(sValues);
    na = numel(alphaValues);

    if abs(constants.shearFactor) < eps
        error('classify_bishop_s_alpha_grid:ZeroShearFactor', ...
            'ShearFactor is zero; provide a nonzero ShearFactor.');
    end
    if abs(constants.alphaFactor) < eps
        error('classify_bishop_s_alpha_grid:ZeroAlphaFactor', ...
            'AlphaFactor is zero; provide a nonzero AlphaFactor.');
    end

    qprimeValues = sValues .* constants.q0 ./ constants.shearFactor;
    pprimeEquationValues = alphaValues ./ constants.alphaFactor;

    if shearMapping == "jperiod"
        shearCoeffs = compute_jperiod_shear_coefficients(geom, constants, opt);
        targetJperiod = constants.JShearFactor .* sValues;
        IprimeValues = NaN(ns, na);
    else
        shearCoeffs = [];
        targetJperiod = NaN(ns, 1);
        IprimeValues = qprimeValues ./ constants.Cq;
    end

    phaseIndices = build_phase_indices(geom, opt);
    lambdaMin = NaN(ns, na);
    phaseIndexMin = NaN(ns, na);
    status = strings(ns, na);
    failures = strings(ns, na);

    if opt.Verbose
        fprintf('classifying s-alpha grid: %d x %d points\n', ns, na);
        fprintf('  q0=%.12e, Cq=%.12e\n', constants.q0, constants.Cq);
        fprintf('  shearFactor=%.12e, alphaFactor=%.12e\n', ...
            constants.shearFactor, constants.alphaFactor);
        fprintf('  shearMapping=%s, JShearFactor=%.12e\n', ...
            shearMapping, constants.JShearFactor);
    end

    for i = 1:ns
        if opt.Verbose
            if shearMapping == "jperiod"
                fprintf('s row %d/%d: s=% .8e, target Jperiod=% .8e\n', ...
                    i, ns, sValues(i), targetJperiod(i));
            else
                fprintf('s row %d/%d: s=% .8e, Iprime=% .8e\n', ...
                    i, ns, sValues(i), IprimeValues(i));
            end
        end

        for j = 1:na
            if shearMapping == "jperiod"
                IprimeCurrent = (targetJperiod(i) ...
                    - shearCoeffs.J0 ...
                    - shearCoeffs.Jp .* pprimeEquationValues(j)) ...
                    ./ shearCoeffs.JIprime;
                IprimeValues(i,j) = IprimeCurrent;
            else
                IprimeCurrent = IprimeValues(i);
            end

            try
                [lambdaMin(i,j), phaseIndexMin(i,j)] = evaluate_phase_scan(geom, ...
                    phaseIndices, constants, opt, IprimeCurrent, pprimeEquationValues(j));
                status(i,j) = classify_lambda(lambdaMin(i,j), opt.StabilityTolerance);
            catch ME
                failures(i,j) = string(ME.message);
                status(i,j) = "failed";
            end

            if opt.Verbose
                if shearMapping == "jperiod"
                    fprintf(['  alpha=% .8e, pprimeEq=% .8e, Iprime=% .8e, ', ...
                        'lambda=% .8e, l0Index=%g, %s\n'], ...
                        alphaValues(j), pprimeEquationValues(j), IprimeCurrent, ...
                        lambdaMin(i,j), phaseIndexMin(i,j), status(i,j));
                else
                    fprintf(['  alpha=% .8e, pprimeEq=% .8e, lambda=% .8e, ', ...
                        'l0Index=%g, %s\n'], ...
                        alphaValues(j), pprimeEquationValues(j), lambdaMin(i,j), ...
                        phaseIndexMin(i,j), status(i,j));
                end
            end
        end
    end

    grid = struct();
    grid.s = sValues;
    grid.alpha = alphaValues;
    grid.qprime = qprimeValues;
    grid.Iprime = IprimeValues;
    grid.pprimeEquation = pprimeEquationValues;
    grid.targetJperiod = targetJperiod;
    grid.shearMapping = char(shearMapping);
    grid.shearCoefficients = shearCoeffs;
    grid.lambdaMin = lambdaMin;
    grid.phaseIndexMin = phaseIndexMin;
    [grid.phaseL0Min, grid.phaseThetaMin] = phase_index_to_l0(geom, phaseIndexMin);
    grid.ballooningPhaseIndices = phaseIndices;
    grid.status = status;
    grid.stable = status == "stable";
    grid.unstable = status == "unstable";
    grid.marginal = status == "marginal";
    grid.failed = status == "failed";
    grid.failures = failures;
    grid.constants = constants;
    grid.options = opt;

    if opt.Plot
        plot_s_alpha_classification(grid);
    end
end

%======================================================================
function phaseIndices = build_phase_indices(geom, opt)
    n = numel(geom.l);
    if ~isempty(opt.BallooningPhaseIndices)
        phaseIndices = unique(round(opt.BallooningPhaseIndices(:))).';
    elseif opt.ScanBallooningPhase
        nphase = min(round(opt.NumBallooningPhases), n);
        phaseIndices = unique(round(linspace(1, n + 1, nphase + 1)));
        phaseIndices = phaseIndices(1:end-1);
    else
        phaseIndices = round(opt.L0Index);
    end

    phaseIndices = phaseIndices(phaseIndices >= 1 & phaseIndices <= n);
    if isempty(phaseIndices)
        error('classify_bishop_s_alpha_grid:NoPhaseIndices', ...
            'No valid ballooning phase indices were supplied.');
    end
end

%======================================================================
function [lambdaBest, phaseBest] = evaluate_phase_scan(geom, phaseIndices, constants, opt, IprimeCurrent, pprimeEquation)
    lambdaBest = Inf;
    phaseBest = phaseIndices(1);
    lastError = [];

    for kk = 1:numel(phaseIndices)
        try
            ch3 = compute_bishop_ch3_terms(geom, ...
                'L0Index', phaseIndices(kk), ...
                'NPeriodsEachSide', opt.NPeriodsEachSide, ...
                'I0', constants.I0, ...
                'Iprime0', IprimeCurrent, ...
                'PprimeEquation', pprimeEquation, ...
                'Plot', false);

            stab = solve_bishop_ballooning_stability(ch3, ...
                'UseExtended', true, ...
                'NumModes', opt.NumModes, ...
                'StabilityTolerance', opt.StabilityTolerance, ...
                'Plot', false);

            if stab.lambdaMin < lambdaBest
                lambdaBest = stab.lambdaMin;
                phaseBest = phaseIndices(kk);
            end
        catch ME
            lastError = ME;
        end
    end

    if ~isfinite(lambdaBest)
        if isempty(lastError)
            error('classify_bishop_s_alpha_grid:PhaseScanFailed', ...
                'No finite lambda was produced by the ballooning phase scan.');
        else
            rethrow(lastError);
        end
    end
end

%======================================================================
function [phaseL0, phaseTheta] = phase_index_to_l0(geom, phaseIndex)
    phaseL0 = NaN(size(phaseIndex));
    phaseTheta = NaN(size(phaseIndex));

    valid = isfinite(phaseIndex);
    if ~any(valid(:))
        return;
    end

    idx = round(phaseIndex(valid));
    l0 = NaN(size(idx));
    theta0 = NaN(size(idx));
    good = idx >= 1 & idx <= numel(geom.l);

    l0(good) = geom.l(idx(good));
    if isfield(geom, 'L') && isfinite(geom.L) && abs(geom.L) > eps
        theta0(good) = 2*pi .* l0(good) ./ geom.L;
    end

    phaseL0(valid) = l0;
    phaseTheta(valid) = theta0;
end

%======================================================================
function constants = compute_constants(geom, opt)
    if isempty(opt.I0)
        I0 = geom.I0;
    else
        I0 = opt.I0;
    end

    l = geom.l(:);
    R = geom.R(:);
    Bp = abs(geom.Bp(:));
    L = geom.L;

    Cq = periodic_integral(l, 1 ./ (R.^2 .* Bp), L) / (2*pi);
    qGeom = I0 * Cq;

    if isempty(opt.Q0)
        if isfield(geom, 'q0') && isfinite(geom.q0)
            q0 = geom.q0;
        else
            q0 = qGeom;
        end
    else
        q0 = opt.Q0;
    end

    if isempty(opt.ShearFactor)
        if isfield(geom, 'psi_relative') && isfinite(geom.psi_relative)
            shearFactor = 2 * geom.psi_relative;
        else
            shearFactor = 1;
        end
    else
        shearFactor = opt.ShearFactor;
    end

    if isempty(opt.AlphaFactor)
        if isfield(geom, 'alphaFactor0') && isfinite(geom.alphaFactor0)
            alphaFactor = geom.alphaFactor0;
        else
            B2ref = mean(Bp.^2 + (I0./R).^2, 'omitnan');
            gradpsiRef = mean(abs(geom.gradpsi), 'omitnan');
            alphaFactor = -2 * q0.^2 .* geom.R0 .* gradpsiRef ./ B2ref;
        end
    else
        alphaFactor = opt.AlphaFactor;
    end

    constants = struct();
    constants.I0 = I0;
    constants.Cq = Cq;
    constants.qGeom = qGeom;
    constants.q0 = q0;
    constants.shearFactor = shearFactor;
    constants.alphaFactor = alphaFactor;

    if isempty(opt.JShearFactor)
        constants.JShearFactor = calibrate_jperiod_per_s(geom, constants, opt);
    else
        constants.JShearFactor = opt.JShearFactor;
    end
end

%======================================================================
function coeffs = compute_jperiod_shear_coefficients(geom, constants, opt)
    J0 = evaluate_jperiod(geom, constants.I0, 0, 0, opt);
    JI = evaluate_jperiod(geom, constants.I0, 1, 0, opt) - J0;
    Jp = evaluate_jperiod(geom, constants.I0, 0, 1, opt) - J0;

    if abs(JI) < 1e-12
        error('classify_bishop_s_alpha_grid:WeakIprimeControl', ...
            'Jperiod is insensitive to Iprime; cannot enforce shear.');
    end

    coeffs = struct();
    coeffs.J0 = J0;
    coeffs.JIprime = JI;
    coeffs.Jpprime = Jp;
    coeffs.Jp = Jp;
end

%======================================================================
function JShearFactor = calibrate_jperiod_per_s(geom, constants, opt)
    if ~isfield(geom, 'shear0') || ~isfinite(geom.shear0) || abs(geom.shear0) < 1e-12
        warning('classify_bishop_s_alpha_grid:NoBaselineShear', ...
            ['geom.shear0 is unavailable. Using JShearFactor=-2*pi; ', ...
             'consider passing JShearFactor explicitly.']);
        JShearFactor = -2*pi;
        return;
    end

    if isfield(geom, 'Iprime0') && isfinite(geom.Iprime0)
        IprimeBase = geom.Iprime0;
    else
        IprimeBase = 0;
    end

    if isfield(geom, 'mu0_pprime0') && isfinite(geom.mu0_pprime0)
        pprimeBase = geom.mu0_pprime0;
    else
        pprimeBase = 0;
    end

    JBase = evaluate_jperiod(geom, constants.I0, IprimeBase, pprimeBase, opt);
    JShearFactor = JBase ./ geom.shear0;
end

%======================================================================
function Jperiod = evaluate_jperiod(geom, I0, Iprime0, pprimeEquation, opt)
    ch3 = compute_bishop_ch3_terms(geom, ...
        'L0Index', opt.L0Index, ...
        'NPeriodsEachSide', 0, ...
        'I0', I0, ...
        'Iprime0', Iprime0, ...
        'PprimeEquation', pprimeEquation, ...
        'Plot', false);

    Jperiod = ch3.eq29.Jperiod;
end

%======================================================================
function val = periodic_integral(l, f, L)
    l = l(:);
    f = f(:);
    lEnd = [l; L];
    fEnd = [f; f(1)];
    val = sum(0.5 .* (fEnd(1:end-1) + fEnd(2:end)) .* diff(lEnd));
end

%======================================================================
function status = classify_lambda(lambda, tol)
    if ~isfinite(lambda)
        status = "failed";
    elseif lambda < -tol
        status = "unstable";
    elseif lambda > tol
        status = "stable";
    else
        status = "marginal";
    end
end

%======================================================================
function plot_s_alpha_classification(grid)
    code = NaN(size(grid.lambdaMin));
    code(grid.stable) = 1;
    code(grid.marginal) = 0;
    code(grid.unstable) = -1;

    figure('Color', 'w', 'Name', 'Bishop s-alpha classification');
    imagesc(grid.alpha, grid.s, code);
    set(gca, 'YDir', 'normal');
    colormap([0.8 0.2 0.2; 0.95 0.95 0.95; 0.2 0.45 0.85]);
    caxis([-1 1]);
    colorbar('Ticks', [-1 0 1], ...
        'TickLabels', {'unstable', 'marginal', 'stable'});
    xlabel('alpha');
    ylabel('s');
    title('Stable/unstable classification');
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
end
