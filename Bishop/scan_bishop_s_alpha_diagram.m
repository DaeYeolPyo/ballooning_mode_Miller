function sa = scan_bishop_s_alpha_diagram(geom, varargin)
%SCAN_BISHOP_S_ALPHA_DIAGRAM Build a Bishop-style s-alpha stability diagram.
%
%   sa = SCAN_BISHOP_S_ALPHA_DIAGRAM(geom)
%   sa = SCAN_BISHOP_S_ALPHA_DIAGRAM(geom, 'Name', value, ...)
%
% This is a marginal-curve wrapper around CLASSIFY_BISHOP_S_ALPHA_GRID.
% The important point is that the requested magnetic shear is held fixed
% while alpha is scanned. By default this is done by solving for I' at each
% (s,alpha) so that the Bishop Eq. (29) secular increment Jperiod matches
% the requested shear.
%
% Definitions used by default
%   s = 2*(psi_GEQDSK - psi_axis)/q * dq/dpsi_GEQDSK
%
%   alpha = AlphaFactor * pprimeEquation
%   pprimeEquation = mu0 * dp/dpsi_GEQDSK
%   AlphaFactor = -2*q^2*R0*<|grad psi|>/<B^2>
%
% Thus a standard GEQDSK equilibrium with pressure decreasing outward has
% pprimeEquation < 0 and alpha > 0.
%
% Stability is evaluated as min_l0 lambda_min(s, alpha, l0). The l0 value
% stored in the output is the Bishop ballooning origin that gives the most
% unstable mode at the corresponding scan or marginal point.

    p = inputParser;
    addParameter(p, 'I0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'IprimeValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'QprimeValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'ShearValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'ShearFactor', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Q0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'AlphaFactor', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'AlphaValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'PprimeScaleValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'PprimeEquationValues', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'ShearMapping', 'jperiod', @(x)ischar(x) || isstring(x));
    addParameter(p, 'JShearFactor', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'L0Index', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'ScanBallooningPhase', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'BallooningPhaseIndices', [], @(x)isnumeric(x) && isvector(x) || isempty(x));
    addParameter(p, 'NumBallooningPhases', 16, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NPeriodsEachSide', 5, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NumModes', 3, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'StabilityTolerance', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'RefineCritical', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'MaxRefineIter', 30, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'ParamTolerance', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Plot', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Verbose', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    constants = compute_local_constants(geom, opt);
    shearScan = build_shear_scan(opt, constants);
    alphaValues = build_alpha_scan(geom, opt, constants);

    if opt.Verbose
        fprintf('s-alpha scan: %d shear points x %d alpha points\n', ...
            numel(shearScan.s), numel(alphaValues));
        fprintf('  q0=%.12e, Cq=%.12e, shearFactor=%.12e\n', ...
            constants.q0, constants.Cq, constants.shearFactor);
        fprintf('  alphaFactor=%.12e, shearMapping=%s\n', ...
            constants.alphaFactor, lower(string(opt.ShearMapping)));
    end

    cls = classify_bishop_s_alpha_grid(geom, ...
        'SValues', shearScan.s, ...
        'AlphaValues', alphaValues, ...
        'I0', opt.I0, ...
        'Q0', opt.Q0, ...
        'ShearFactor', opt.ShearFactor, ...
        'AlphaFactor', opt.AlphaFactor, ...
        'ShearMapping', opt.ShearMapping, ...
        'JShearFactor', opt.JShearFactor, ...
        'L0Index', opt.L0Index, ...
        'ScanBallooningPhase', opt.ScanBallooningPhase, ...
        'BallooningPhaseIndices', opt.BallooningPhaseIndices, ...
        'NumBallooningPhases', opt.NumBallooningPhases, ...
        'NPeriodsEachSide', opt.NPeriodsEachSide, ...
        'NumModes', opt.NumModes, ...
        'StabilityTolerance', opt.StabilityTolerance, ...
        'Plot', false, ...
        'Verbose', opt.Verbose);

    n = numel(cls.s);
    scans = cell(n, 1);
    alphaRoots = cell(n, 1);
    pprimeRoots = cell(n, 1);
    alphaRootPhaseIndex = cell(n, 1);
    alphaRootL0 = cell(n, 1);
    alphaRootTheta = cell(n, 1);
    firstAlpha = NaN(n, 1);
    firstPprimeEquation = NaN(n, 1);
    firstPhaseIndex = NaN(n, 1);
    firstL0 = NaN(n, 1);
    firstTheta = NaN(n, 1);
    nRoots = zeros(n, 1);

    for i = 1:n
        scans{i} = struct( ...
            'alpha', cls.alpha(:), ...
            'pprimeEquation', cls.pprimeEquation(:), ...
            'Iprime', cls.Iprime(i,:).', ...
            'lambdaMin', cls.lambdaMin(i,:).', ...
            'phaseIndex', cls.phaseIndexMin(i,:).', ...
            'phaseL0', cls.phaseL0Min(i,:).', ...
            'phaseTheta', cls.phaseThetaMin(i,:).', ...
            'status', cls.status(i,:).');

        roots = find_alpha_roots(geom, cls.s(i), cls.alpha(:), ...
            cls.lambdaMin(i,:).', cls.phaseIndexMin(i,:).', ...
            cls.phaseL0Min(i,:).', cls.phaseThetaMin(i,:).', ...
            cls.constants, opt);

        nRoots(i) = numel(roots);
        if isempty(roots)
            alphaRoots{i} = [];
            pprimeRoots{i} = [];
            alphaRootPhaseIndex{i} = [];
            alphaRootL0{i} = [];
            alphaRootTheta{i} = [];
        else
            aa = arrayfun(@(r) r.alpha, roots(:));
            pp = aa(:) ./ cls.constants.alphaFactor;
            phaseIdx = arrayfun(@(r) r.phaseIndex, roots(:));
            phaseL0 = arrayfun(@(r) r.phaseL0, roots(:));
            phaseTheta = arrayfun(@(r) r.phaseTheta, roots(:));
            alphaRoots{i} = aa(:);
            pprimeRoots{i} = pp(:);
            alphaRootPhaseIndex{i} = phaseIdx(:);
            alphaRootL0{i} = phaseL0(:);
            alphaRootTheta{i} = phaseTheta(:);
            firstAlpha(i) = aa(1);
            firstPprimeEquation(i) = pp(1);
            firstPhaseIndex(i) = phaseIdx(1);
            firstL0(i) = phaseL0(1);
            firstTheta(i) = phaseTheta(1);
        end
    end

    sa = struct();
    sa.s = cls.s(:);
    sa.alpha = cls.alpha(:).';
    sa.Iprime = cls.Iprime;
    sa.qprime = cls.qprime(:);
    sa.q0 = cls.constants.q0;
    sa.Cq = cls.constants.Cq;
    sa.shearFactor = cls.constants.shearFactor;
    sa.alphaFactor = cls.constants.alphaFactor;
    sa.shearMapping = cls.shearMapping;
    sa.targetJperiod = cls.targetJperiod;
    sa.alphaRoots = alphaRoots;
    sa.pprimeEquationRoots = pprimeRoots;
    sa.alphaRootPhaseIndex = alphaRootPhaseIndex;
    sa.alphaRootL0 = alphaRootL0;
    sa.alphaRootTheta = alphaRootTheta;
    sa.firstAlpha = firstAlpha;
    sa.firstPprimeEquation = firstPprimeEquation;
    sa.firstPhaseIndex = firstPhaseIndex;
    sa.firstL0 = firstL0;
    sa.firstTheta = firstTheta;
    sa.nRoots = nRoots;
    sa.phaseIndexMin = cls.phaseIndexMin;
    sa.phaseL0Min = cls.phaseL0Min;
    sa.phaseThetaMin = cls.phaseThetaMin;
    sa.ballooningPhaseIndices = cls.ballooningPhaseIndices;
    sa.scans = scans;
    sa.classification = cls;
    sa.options = opt;

    if opt.Plot
        plot_s_alpha_scan(sa);
    end
end

%======================================================================
function constants = compute_local_constants(geom, opt)
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
end

%======================================================================
function shearScan = build_shear_scan(opt, constants)
    supplied = [~isempty(opt.IprimeValues), ~isempty(opt.QprimeValues), ~isempty(opt.ShearValues)];
    if sum(supplied) > 1
        error('scan_bishop_s_alpha_diagram:TooManyShearInputs', ...
            'Supply only one of IprimeValues, QprimeValues, or ShearValues.');
    end

    if ~isempty(opt.IprimeValues)
        IprimeLabel = opt.IprimeValues(:);
        qprime = constants.Cq .* IprimeLabel;
        s = constants.shearFactor .* qprime ./ constants.q0;
    elseif ~isempty(opt.QprimeValues)
        qprime = opt.QprimeValues(:);
        IprimeLabel = qprime ./ constants.Cq;
        s = constants.shearFactor .* qprime ./ constants.q0;
    elseif ~isempty(opt.ShearValues)
        s = opt.ShearValues(:);
        qprime = s .* constants.q0 ./ constants.shearFactor;
        IprimeLabel = qprime ./ constants.Cq;
    else
        s = linspace(-2, 2, 9).';
        qprime = s .* constants.q0 ./ constants.shearFactor;
        IprimeLabel = qprime ./ constants.Cq;
    end

    shearScan = struct();
    shearScan.s = s(:);
    shearScan.qprime = qprime(:);
    shearScan.IprimeLabel = IprimeLabel(:);
end

%======================================================================
function alphaValues = build_alpha_scan(geom, opt, constants)
    if ~isempty(opt.AlphaValues)
        alphaValues = opt.AlphaValues(:).';
    elseif ~isempty(opt.PprimeEquationValues)
        alphaValues = constants.alphaFactor .* opt.PprimeEquationValues(:).';
    elseif ~isempty(opt.PprimeScaleValues)
        if isfield(geom, 'mu0_pprime0')
            basePprime = geom.mu0_pprime0;
        elseif isfield(geom, 'pprime0')
            basePprime = 4*pi*1e-7 * geom.pprime0;
        else
            error('scan_bishop_s_alpha_diagram:MissingBasePprime', ...
                'Need geom.mu0_pprime0 or AlphaValues.');
        end
        alphaValues = constants.alphaFactor .* basePprime .* opt.PprimeScaleValues(:).';
    else
        alphaValues = linspace(0, 3, 13);
    end

    alphaValues = unique(alphaValues(isfinite(alphaValues)), 'sorted');
    if numel(alphaValues) < 2
        error('scan_bishop_s_alpha_diagram:TooFewAlphaPoints', ...
            'At least two finite alpha points are required.');
    end
end

%======================================================================
function roots = find_alpha_roots(geom, sValue, alphaValues, lambdaValues, ...
    phaseIndexValues, phaseL0Values, phaseThetaValues, constants, opt)
    roots = struct('alpha', {}, 'lambdaMin', {}, ...
        'phaseIndex', {}, 'phaseL0', {}, 'phaseTheta', {}, ...
        'leftIndex', {}, 'rightIndex', {}, ...
        'leftAlpha', {}, 'rightAlpha', {});

    finite = isfinite(lambdaValues);
    tol = opt.StabilityTolerance;

    for i = 1:(numel(alphaValues)-1)
        if ~finite(i) || ~finite(i+1)
            continue;
        end

        f1 = lambdaValues(i);
        f2 = lambdaValues(i+1);

        if abs(f1) <= tol
            roots(end+1) = make_root(alphaValues(i), f1, ...
                phaseIndexValues(i), phaseL0Values(i), phaseThetaValues(i), i, i, ...
                alphaValues(i), alphaValues(i)); %#ok<AGROW>
            continue;
        end

        if f1 * f2 > 0
            continue;
        end

        if abs(f2) <= tol
            roots(end+1) = make_root(alphaValues(i+1), f2, ...
                phaseIndexValues(i+1), phaseL0Values(i+1), phaseThetaValues(i+1), i+1, i+1, ...
                alphaValues(i+1), alphaValues(i+1)); %#ok<AGROW>
            continue;
        end

        if opt.RefineCritical
            [aRoot, fRoot, phaseRoot] = refine_alpha_root(geom, sValue, ...
                alphaValues(i), alphaValues(i+1), f1, f2, constants, opt);
        else
            aRoot = 0.5 * (alphaValues(i) + alphaValues(i+1));
            [fRoot, phaseRoot] = evaluate_lambda_s_alpha(geom, sValue, aRoot, constants, opt);
        end

        roots(end+1) = make_root(aRoot, fRoot, ...
            phaseRoot.phaseIndex, phaseRoot.phaseL0, phaseRoot.phaseTheta, i, i+1, ...
            alphaValues(i), alphaValues(i+1)); %#ok<AGROW>
    end
end

%======================================================================
function root = make_root(alpha, lambdaMin, phaseIndex, phaseL0, phaseTheta, ...
    leftIndex, rightIndex, leftAlpha, rightAlpha)
    root = struct();
    root.alpha = alpha;
    root.lambdaMin = lambdaMin;
    root.phaseIndex = phaseIndex;
    root.phaseL0 = phaseL0;
    root.phaseTheta = phaseTheta;
    root.leftIndex = leftIndex;
    root.rightIndex = rightIndex;
    root.leftAlpha = leftAlpha;
    root.rightAlpha = rightAlpha;
end

%======================================================================
function [aRoot, fRoot, phaseRoot] = refine_alpha_root(geom, sValue, a1, a2, f1, f2, constants, opt)
    a = a1;
    b = a2;
    fa = f1;

    if f1 * f2 > 0
        error('scan_bishop_s_alpha_diagram:BadBracket', ...
            'Refinement requires a sign-changing bracket.');
    end

    for it = 1:round(opt.MaxRefineIter)
        c = 0.5 * (a + b);
        [fc, phaseC] = evaluate_lambda_s_alpha(geom, sValue, c, constants, opt);

        if abs(fc) <= opt.StabilityTolerance || ...
           abs(b - a) <= opt.ParamTolerance * max(1, max(abs([a, b])))
            aRoot = c;
            fRoot = fc;
            phaseRoot = phaseC;
            return;
        end

        if fa * fc <= 0
            b = c;
        else
            a = c;
            fa = fc;
        end
    end

    aRoot = 0.5 * (a + b);
    [fRoot, phaseRoot] = evaluate_lambda_s_alpha(geom, sValue, aRoot, constants, opt);
end

%======================================================================
function [lambdaMin, phaseInfo] = evaluate_lambda_s_alpha(geom, sValue, alphaValue, constants, opt)
    cls = classify_bishop_s_alpha_grid(geom, ...
        'SValues', sValue, ...
        'AlphaValues', alphaValue, ...
        'I0', opt.I0, ...
        'Q0', constants.q0, ...
        'ShearFactor', constants.shearFactor, ...
        'AlphaFactor', constants.alphaFactor, ...
        'ShearMapping', opt.ShearMapping, ...
        'JShearFactor', opt.JShearFactor, ...
        'L0Index', opt.L0Index, ...
        'ScanBallooningPhase', opt.ScanBallooningPhase, ...
        'BallooningPhaseIndices', opt.BallooningPhaseIndices, ...
        'NumBallooningPhases', opt.NumBallooningPhases, ...
        'NPeriodsEachSide', opt.NPeriodsEachSide, ...
        'NumModes', opt.NumModes, ...
        'StabilityTolerance', opt.StabilityTolerance, ...
        'Plot', false, ...
        'Verbose', false);

    lambdaMin = cls.lambdaMin(1,1);
    phaseInfo = make_phase_info(cls.phaseIndexMin(1,1), ...
        cls.phaseL0Min(1,1), cls.phaseThetaMin(1,1));
end

%======================================================================
function phaseInfo = make_phase_info(phaseIndex, phaseL0, phaseTheta)
    phaseInfo = struct();
    phaseInfo.phaseIndex = phaseIndex;
    phaseInfo.phaseL0 = phaseL0;
    phaseInfo.phaseTheta = phaseTheta;
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
function plot_s_alpha_scan(sa)
    figure('Color', 'w', 'Name', 'Bishop s-alpha diagram');
    hold on;

    imagesc(sa.alpha, sa.s, double(sa.classification.unstable));
    set(gca, 'YDir', 'normal');
    colormap([0.2 0.45 0.85; 0.8 0.2 0.2]);
    caxis([0 1]);
    colorbar('Ticks', [0 1], 'TickLabels', {'stable', 'unstable'});

    for i = 1:numel(sa.s)
        aa = sa.alphaRoots{i};
        if isempty(aa)
            continue;
        end
        plot(aa, sa.s(i) .* ones(size(aa)), 'ko', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 5);
    end

    good = isfinite(sa.firstAlpha);
    if any(good)
        plot(sa.firstAlpha(good), sa.s(good), 'k-', 'LineWidth', 1.5);
    end

    xlabel('alpha');
    ylabel('s');
    title('Bishop s-alpha stable/unstable scan');
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
end
