function scan = scan_bishop_critical_pprime(geom, varargin)
%SCAN_BISHOP_CRITICAL_PPRIME Scan p' and find marginal ballooning stability.
%
%   scan = SCAN_BISHOP_CRITICAL_PPRIME(geom)
%   scan = SCAN_BISHOP_CRITICAL_PPRIME(geom, 'Name', value, ...)
%
% This treats the pressure-gradient value used in Bishop Chapter 3 as a
% scalar parameter and solves
%
%   lambda_min(pprimeEquation) = 0
%
% for the marginal/critical pressure gradient.
%
% Input
%   geom : output of build_bishop_geometry(...)
%
% Name-value options
%   'BasePprimeEquation'   : baseline value. Default geom.mu0_pprime0.
%                            With build_bishop_geometry this is
%                            mu0*dp/dpsi_GEQDSK.
%   'ScaleValues'          : values s used in pprimeEquation = s*base.
%                            Default linspace(0, 2, 17).
%   'PprimeEquationValues' : direct values of pprimeEquation. If supplied,
%                            ScaleValues is ignored.
%   'L0Index'              : reference index for Bishop integral when
%                            ScanBallooningPhase=false. Default 1.
%   'ScanBallooningPhase'  : true scans l0/L0Index and keeps the most
%                            unstable eigenvalue. Default true.
%   'BallooningPhaseIndices': explicit L0Index values for the l0 scan.
%   'NumBallooningPhases'  : number of approximately uniform phases if
%                            ScanBallooningPhase=true and explicit indices
%                            are not supplied. Default 16.
%   'I0', 'Iprime0'        : optional overrides passed to
%                            compute_bishop_ch3_terms.
%   'NPeriodsEachSide'     : ballooning theta half-window. Default 5.
%   'NumModes'             : eigenmodes requested from stability solve.
%                            Default 3.
%   'StabilityTolerance'   : lambda tolerance near zero. Default 1e-8.
%   'RefineCritical'       : refine sign-change brackets by bisection.
%                            Default true.
%   'MaxRefineIter'        : bisection iterations. Default 35.
%   'ParamTolerance'       : bisection parameter tolerance. Default 1e-8.
%   'Plot'                 : true/false. Default true.
%   'Verbose'              : true/false. Default true.
%
% Output
%   scan.parameter         : scale values or direct pprimeEquation values
%   scan.pprimeEquation    : values passed to compute_bishop_ch3_terms
%   scan.lambdaMin         : lowest eigenvalue at each scan point
%   scan.status            : stable/marginal/unstable at each scan point
%   scan.roots             : all refined sign-change roots
%   scan.critical          : first refined root in scan order, if found
%
% Notes
%   pprimeEquation is mu0*dp/dpsi_GEQDSK by default when geom comes from
%   build_bishop_geometry. A standard pressure-decreasing equilibrium has
%   pprimeEquation < 0.

    p = inputParser;
    addParameter(p, 'BasePprimeEquation', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'ScaleValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'PprimeEquationValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'I0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Iprime0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'L0Index', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'ScanBallooningPhase', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'BallooningPhaseIndices', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'NumBallooningPhases', 16, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NPeriodsEachSide', 5, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NumModes', 3, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'StabilityTolerance', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'RefineCritical', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'MaxRefineIter', 35, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'ParamTolerance', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Plot', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Verbose', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    [parameter, pprimeEquation, parameterName, basePprimeEquation] = ...
        build_scan_grid(geom, opt);

    nscan = numel(parameter);
    lambdaMin = NaN(nscan, 1);
    phaseIndexMin = NaN(nscan, 1);
    phaseL0Min = NaN(nscan, 1);
    phaseThetaMin = NaN(nscan, 1);
    status = strings(nscan, 1);
    failures = strings(nscan, 1);

    if opt.Verbose
        fprintf('pprime scan: %d points, %s range [%.8g, %.8g]\n', ...
            nscan, parameterName, min(parameter), max(parameter));
        fprintf('base pprimeEquation = %.12e\n', basePprimeEquation);
    end

    for i = 1:nscan
        try
            [lambdaMin(i), phaseInfo] = evaluate_lambda(geom, pprimeEquation(i), opt);
            phaseIndexMin(i) = phaseInfo.phaseIndex;
            phaseL0Min(i) = phaseInfo.phaseL0;
            phaseThetaMin(i) = phaseInfo.phaseTheta;
            status(i) = classify_lambda(lambdaMin(i), opt.StabilityTolerance);
        catch ME
            failures(i) = string(ME.message);
            status(i) = "failed";
        end

        if opt.Verbose
            fprintf(['  %s=% .8e, pprimeEq=% .8e, lambda=% .12e, ', ...
                'l0Index=%g, %s\n'], ...
                parameterName, parameter(i), pprimeEquation(i), lambdaMin(i), ...
                phaseIndexMin(i), status(i));
        end
    end

    roots = find_and_refine_roots(geom, parameter, pprimeEquation, ...
        lambdaMin, phaseIndexMin, phaseL0Min, phaseThetaMin, parameterName, opt);

    if isempty(roots)
        critical = [];
        if opt.Verbose
            fprintf('  no marginal crossing found in this scan range\n');
        end
    else
        critical = roots(1);
        if opt.Verbose
            for i = 1:numel(roots)
                fprintf(['  root %d: %s=%.12e, pprimeEq=%.12e, ', ...
                    'lambda=%.12e, l0Index=%g\n'], ...
                    i, parameterName, roots(i).parameter, ...
                    roots(i).pprimeEquation, roots(i).lambdaMin, ...
                    roots(i).phaseIndex);
            end
        end
    end

    scan = struct();
    scan.parameterName = parameterName;
    scan.parameter = parameter(:);
    scan.pprimeEquation = pprimeEquation(:);
    scan.basePprimeEquation = basePprimeEquation;
    scan.lambdaMin = lambdaMin;
    scan.phaseIndexMin = phaseIndexMin;
    scan.phaseL0Min = phaseL0Min;
    scan.phaseThetaMin = phaseThetaMin;
    scan.ballooningPhaseIndices = build_phase_indices(geom, opt);
    scan.status = status;
    scan.failures = failures;
    scan.roots = roots;
    scan.critical = critical;
    scan.options = opt;

    if opt.Plot
        plot_pprime_scan(scan);
    end
end

%======================================================================
function [parameter, pprimeEquation, parameterName, basePprimeEquation] = ...
    build_scan_grid(geom, opt)

    if isempty(opt.BasePprimeEquation)
        if isfield(geom, 'mu0_pprime0')
            basePprimeEquation = geom.mu0_pprime0;
        elseif isfield(geom, 'pprime0')
            basePprimeEquation = 4*pi*1e-7 * geom.pprime0;
        else
            error('scan_bishop_critical_pprime:MissingBasePprime', ...
                ['BasePprimeEquation was not supplied and geom has no ', ...
                 'mu0_pprime0/pprime0 field.']);
        end
    else
        basePprimeEquation = opt.BasePprimeEquation;
    end

    if ~isempty(opt.PprimeEquationValues)
        pprimeEquation = opt.PprimeEquationValues(:);
        parameter = pprimeEquation;
        parameterName = 'pprimeEquation';
    else
        if abs(basePprimeEquation) < eps
            error('scan_bishop_critical_pprime:ZeroBasePprime', ...
                ['BasePprimeEquation is zero. Supply PprimeEquationValues ', ...
                 'for a direct scan.']);
        end

        if isempty(opt.ScaleValues)
            scaleValues = linspace(0, 2, 17).';
        else
            scaleValues = opt.ScaleValues(:);
        end

        parameter = scaleValues;
        pprimeEquation = basePprimeEquation .* scaleValues;
        parameterName = 'scale';
    end

    if numel(parameter) < 2
        error('scan_bishop_critical_pprime:TooFewScanPoints', ...
            'At least two scan points are required.');
    end
end

%======================================================================
function [lambdaMin, phaseInfo] = evaluate_lambda(geom, pprimeEquation, opt)
    phaseIndices = build_phase_indices(geom, opt);
    [lambdaMin, phaseIndex] = evaluate_phase_scan(geom, phaseIndices, pprimeEquation, opt);
    [phaseL0, phaseTheta] = phase_index_to_l0(geom, phaseIndex);
    phaseInfo = make_phase_info(phaseIndex, phaseL0, phaseTheta);
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
        error('scan_bishop_critical_pprime:NoPhaseIndices', ...
            'No valid ballooning phase indices were supplied.');
    end
end

%======================================================================
function [lambdaBest, phaseBest] = evaluate_phase_scan(geom, phaseIndices, pprimeEquation, opt)
    lambdaBest = Inf;
    phaseBest = phaseIndices(1);
    lastError = [];

    for kk = 1:numel(phaseIndices)
        try
            ch3 = compute_bishop_ch3_terms(geom, ...
                'L0Index', phaseIndices(kk), ...
                'NPeriodsEachSide', opt.NPeriodsEachSide, ...
                'PprimeEquation', pprimeEquation, ...
                'I0', opt.I0, ...
                'Iprime0', opt.Iprime0, ...
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
            error('scan_bishop_critical_pprime:PhaseScanFailed', ...
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
function phaseInfo = make_phase_info(phaseIndex, phaseL0, phaseTheta)
    phaseInfo = struct();
    phaseInfo.phaseIndex = phaseIndex;
    phaseInfo.phaseL0 = phaseL0;
    phaseInfo.phaseTheta = phaseTheta;
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
function roots = find_and_refine_roots(geom, parameter, pprimeEquation, ...
    lambdaMin, phaseIndexMin, phaseL0Min, phaseThetaMin, parameterName, opt)

    roots = struct( ...
        'parameterName', {}, ...
        'parameter', {}, ...
        'pprimeEquation', {}, ...
        'lambdaMin', {}, ...
        'phaseIndex', {}, ...
        'phaseL0', {}, ...
        'phaseTheta', {}, ...
        'leftIndex', {}, ...
        'rightIndex', {}, ...
        'leftParameter', {}, ...
        'rightParameter', {});

    finite = isfinite(lambdaMin);
    tol = opt.StabilityTolerance;

    for i = 1:(numel(parameter)-1)
        if ~finite(i) || ~finite(i+1)
            continue;
        end

        f1 = lambdaMin(i);
        f2 = lambdaMin(i+1);

        if abs(f1) <= tol
            root = make_root(parameterName, parameter(i), pprimeEquation(i), ...
                f1, phaseIndexMin(i), phaseL0Min(i), phaseThetaMin(i), ...
                i, i, parameter(i), parameter(i));
            roots(end+1) = root; %#ok<AGROW>
            continue;
        end

        if f1 * f2 > 0
            continue;
        end

        if abs(f2) <= tol
            root = make_root(parameterName, parameter(i+1), pprimeEquation(i+1), ...
                f2, phaseIndexMin(i+1), phaseL0Min(i+1), phaseThetaMin(i+1), ...
                i+1, i+1, parameter(i+1), parameter(i+1));
            roots(end+1) = root; %#ok<AGROW>
            continue;
        end

        if opt.RefineCritical
            [pRoot, ppRoot, fRoot, phaseRoot] = refine_bracket(geom, ...
                parameter(i), parameter(i+1), ...
                pprimeEquation(i), pprimeEquation(i+1), ...
                f1, f2, opt);
        else
            pRoot = 0.5 * (parameter(i) + parameter(i+1));
            ppRoot = 0.5 * (pprimeEquation(i) + pprimeEquation(i+1));
            [fRoot, phaseRoot] = evaluate_lambda(geom, ppRoot, opt);
        end

        root = make_root(parameterName, pRoot, ppRoot, fRoot, ...
            phaseRoot.phaseIndex, phaseRoot.phaseL0, phaseRoot.phaseTheta, ...
            i, i+1, parameter(i), parameter(i+1));
        roots(end+1) = root; %#ok<AGROW>
    end
end

%======================================================================
function root = make_root(parameterName, parameter, pprimeEquation, ...
    lambdaMin, phaseIndex, phaseL0, phaseTheta, ...
    leftIndex, rightIndex, leftParameter, rightParameter)

    root = struct();
    root.parameterName = parameterName;
    root.parameter = parameter;
    root.pprimeEquation = pprimeEquation;
    root.lambdaMin = lambdaMin;
    root.phaseIndex = phaseIndex;
    root.phaseL0 = phaseL0;
    root.phaseTheta = phaseTheta;
    root.leftIndex = leftIndex;
    root.rightIndex = rightIndex;
    root.leftParameter = leftParameter;
    root.rightParameter = rightParameter;
end

%======================================================================
function [pRoot, ppRoot, fRoot, phaseRoot] = refine_bracket(geom, p1, p2, pp1, pp2, f1, f2, opt)
    a = p1;
    b = p2;
    fa = f1;
    fb = f2;

    if fa * fb > 0
        error('scan_bishop_critical_pprime:BadBracket', ...
            'Refinement requires a sign-changing bracket.');
    end

    for it = 1:round(opt.MaxRefineIter)
        c = 0.5 * (a + b);
        ppC = pp1 + (pp2 - pp1) .* (c - p1) ./ (p2 - p1);
        [fc, phaseC] = evaluate_lambda(geom, ppC, opt);

        if abs(fc) <= opt.StabilityTolerance || ...
           abs(b - a) <= opt.ParamTolerance * max(1, max(abs([a, b])))
            pRoot = c;
            ppRoot = ppC;
            fRoot = fc;
            phaseRoot = phaseC;
            return;
        end

        if fa * fc <= 0
            b = c;
            fb = fc; %#ok<NASGU>
        else
            a = c;
            fa = fc;
        end
    end

    pRoot = 0.5 * (a + b);
    ppRoot = pp1 + (pp2 - pp1) .* (pRoot - p1) ./ (p2 - p1);
    [fRoot, phaseRoot] = evaluate_lambda(geom, ppRoot, opt);
end

%======================================================================
function plot_pprime_scan(scan)
    figure('Color', 'w', 'Name', 'Bishop pprime critical scan');

    plot(scan.parameter, scan.lambdaMin, 'o-', 'LineWidth', 1.4);
    hold on;
    yline(0, 'k--');
    grid on;
    xlabel(scan.parameterName);
    ylabel('lambda_min');
    title('Critical pressure-gradient scan');

    if ~isempty(scan.critical)
        plot(scan.critical.parameter, scan.critical.lambdaMin, ...
            'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
        legend('scan', 'marginal', 'critical', 'Location', 'best');
    else
        legend('scan', 'marginal', 'Location', 'best');
    end
end
