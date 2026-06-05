function ch3 = compute_bishop_ch3_terms(geom, varargin)
%COMPUTE_BISHOP_CH3_TERMS Compute Bishop Chapter 3 ballooning terms.
%
%   ch3 = COMPUTE_BISHOP_CH3_TERMS(geom)
%   ch3 = COMPUTE_BISHOP_CH3_TERMS(geom, 'Name', value, ...)
%
% Input
%   geom : output of build_bishop_geometry(...)
%
% Name-value options
%   'L0Index'           : reference index for the ballooning integral.
%                         Default 1.
%   'NPeriodsEachSide'  : if >0, also return ch3.extended on
%                         theta in approximately
%                         [-2*pi*N, 2*pi*N). Default 0.
%   'UseMu0Pprime'      : true uses geom.mu0_pprime0 = mu0*dp/dpsi_GEQDSK.
%                         Default true.
%   'PprimeEquation'    : explicit pressure-gradient value used in the
%                         Bishop equations. Default [].
%   'DriveFactor'       : multiplier for Eq. (31). Bishop Eq. (31) is
%                         2*p'*(B x grad S).kappa/B^2; the compact bracket
%                         below is multiplied by this factor. Default 2.
%   'I0', 'Iprime0'     : optional overrides for I and dI/dpsi. Default [].
%   'MinBp'             : lower bound for Bp. Default 1e-12.
%   'MinCurvatureRadius': lower bound for |Rc|. Default 1e-10.
%   'ModeF'             : optional ballooning amplitude F(theta). If
%                         supplied, equation terms are evaluated.
%   'Plot'              : true/false. Default false.
%
% Output structure
%   ch3.eq29 : pieces of Bishop Eq. (29), |grad S|^2.
%   ch3.eq31 : pieces of Bishop Eq. (31), curvature-pressure drive.
%   ch3.eq32 : coefficients in
%
%       Bp d/dl { A_inner dF/dl } + C_drive F = 0,
%
%       A_inner = (|grad S|^2/B^2) Bp.
%
%   For a standard self-adjoint form divided by Bp:
%
%       d/dl { A_inner dF/dl } + C_drive/Bp F = 0.
%
% Ballooning extension
%   Periodic geometric quantities are repeated in theta, but the integral
%   J(theta) in Eq. (29)/(31) is continued as
%
%       J(theta + 2*pi*k) = J(theta) + k*J_period.
%
%   This is the part that carries the extended ballooning representation.

    p = inputParser;
    addParameter(p, 'L0Index', 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'NPeriodsEachSide', 0, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'UseMu0Pprime', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'PprimeEquation', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'I0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Iprime0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'DriveFactor', 2, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'MinBp', 1e-12, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MinCurvatureRadius', 1e-10, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'ModeF', [], @(x)isnumeric(x) || isempty(x));
    addParameter(p, 'Plot', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    data = prepare_one_period_data(geom, opt);
    one = compute_terms_on_path(data, opt);

    ch3 = one;
    ch3.options = opt;

    nside = round(opt.NPeriodsEachSide);
    if nside > 0
        ch3.extended = extend_ch3_terms(one, nside, opt);
    else
        ch3.extended = [];
    end

    if ~isempty(opt.ModeF)
        if nside > 0
            ch3.modeTerms = evaluate_mode_terms(ch3.extended, opt.ModeF(:));
        else
            ch3.modeTerms = evaluate_mode_terms(ch3, opt.ModeF(:));
        end
    else
        ch3.modeTerms = [];
    end

    if opt.Plot
        plot_ch3_terms(ch3);
    end
end

%======================================================================
function data = prepare_one_period_data(geom, opt)
    required = {'l','L','R','Z','u','h0','Bp','R0','I0','Iprime0'};
    for i = 1:numel(required)
        if ~isfield(geom, required{i})
            error('compute_bishop_ch3_terms:MissingField', ...
                'geom.%s is required.', required{i});
        end
    end

    l = geom.l(:);
    n = numel(l);
    i0 = round(opt.L0Index);
    if i0 < 1 || i0 > n
        error('compute_bishop_ch3_terms:BadL0Index', ...
            'L0Index must be between 1 and %d.', n);
    end

    idx = [i0:n, 1:i0-1];
    l0 = l(i0);
    lShift = l(idx) - l0;
    lShift(lShift < 0) = lShift(lShift < 0) + geom.L;

    R = geom.R(idx);
    Z = geom.Z(idx);
    u = unwrap(geom.u(idx));
    Bp = max(abs(geom.Bp(idx)), opt.MinBp);

    if isfield(geom, 'Rc_signed')
        Rc = geom.Rc_signed(idx);
    elseif isfield(geom, 'Rc')
        Rc = geom.Rc(idx);
    else
        Rc = 1 ./ geom.kappa_signed(idx);
    end

    smallRc = abs(Rc) < opt.MinCurvatureRadius | ~isfinite(Rc);
    if any(smallRc)
        s = sign(Rc(smallRc));
        s(s == 0) = 1;
        Rc(smallRc) = s .* opt.MinCurvatureRadius;
    end

    if isempty(opt.I0)
        I0 = geom.I0;
    else
        I0 = opt.I0;
    end

    if isempty(opt.Iprime0)
        Iprime0 = geom.Iprime0;
    else
        Iprime0 = opt.Iprime0;
    end

    if isempty(opt.PprimeEquation)
        if opt.UseMu0Pprime
            if isfield(geom, 'mu0_pprime0')
                pprimeEq = geom.mu0_pprime0;
            else
                pprimeEq = 4*pi*1e-7 * geom.pprime0;
            end
        else
            pprimeEq = geom.pprime0;
        end
    else
        pprimeEq = opt.PprimeEquation;
    end

    if isfield(geom, 'x0') && isfinite(geom.x0)
        x0 = geom.x0;
    else
        x0 = geom.R(1);
    end
    h0 = R ./ x0;

    data = struct();
    data.l = lShift(:);
    data.L = geom.L;
    data.theta = 2*pi * lShift(:) ./ geom.L;
    data.R = R(:);
    data.Z = Z(:);
    data.u = u(:);
    data.h0 = h0(:);
    data.x0 = x0;
    data.Bp = Bp(:);
    data.Rc = Rc(:);
    data.I0 = I0;
    data.Iprime0 = Iprime0;
    data.pprimeEquation = pprimeEq;

    if isfield(geom, 'pprime0')
        data.pprime0 = geom.pprime0;
    else
        data.pprime0 = NaN;
    end

    if isfield(geom, 'mu0_pprime0')
        data.mu0_pprime0 = geom.mu0_pprime0;
    else
        data.mu0_pprime0 = NaN;
    end

    if isfield(geom, 'direction')
        data.direction = geom.direction;
    else
        data.direction = 'clockwise';
    end
end

%======================================================================
function out = compute_terms_on_path(data, opt)
    l = data.l;
    L = data.L;
    R = data.R;
    u = data.u;
    h0 = data.h0;
    x0 = data.x0;
    Bp = data.Bp;
    Rc = data.Rc;
    I0 = data.I0;
    Iprime0 = data.Iprime0;
    pprimeEq = data.pprimeEquation;

    Btor = I0 ./ R;
    B2 = Bp.^2 + Btor.^2;
    B4 = B2.^2;
    dB2dl = periodic_derivative(l, B2, L);

    metricFactor = 1 + I0^2 ./ (x0^2 .* h0.^2 .* Bp.^2);

    % Bishop Eq. (29) integral. With build_bishop_geometry defaults,
    % pprimeEq is mu0*dp/dpsi_GEQDSK.
    Q29_Iprime = (Iprime0 / I0) .* metricFactor;
    Q29_pressure = pprimeEq ./ Bp.^2;
    Q29_curvature = (2 ./ Rc) .* ...
        (1 + Rc .* sin(u) ./ (x0 .* h0));
    Q29_poloidalField = -1 ./ (x0 .* h0 .* Bp);

    Q29 = Q29_Iprime + Q29_pressure + Q29_curvature + Q29_poloidalField;
    integrandJ = Q29 ./ (h0.^2 .* Bp);
    [J, Jperiod] = periodic_cumtrapz(l, integrandJ, L);

    gradS2_metric = metricFactor ./ (x0^2 .* h0.^2);
    gradS2_shear = (h0.^2 .* Bp.^2 .* I0^2 ./ x0^2) .* J.^2;
    gradS2 = gradS2_metric + gradS2_shear;

    % Bishop Eq. (31) bracket and pressure-curvature drive.
    curv_toroidicity = 2 .* I0^2 .* sin(u) ./ ...
        (x0^3 .* h0.^3 .* B2 .* Bp);
    curv_normal = -2 .* Bp ./ (Rc .* B2);
    curv_shear = -I0^2 .* h0 .* Bp ./ (x0 .* B4) .* dB2dl .* J;
    curvBracket = curv_toroidicity + curv_normal + curv_shear;
    curvDrive = opt.DriveFactor .* pprimeEq ./ (x0 .* h0) .* curvBracket;

    A_inner = gradS2 ./ B2 .* Bp;
    C_drive = curvDrive;
    C_dividedByBp = C_drive ./ Bp;

    out = data;

    out.Btor = Btor;
    out.B2 = B2;
    out.B4 = B4;
    out.dB2dl = dB2dl;

    out.eq29 = struct();
    out.eq29.metricFactor = metricFactor;
    out.eq29.Q_Iprime = Q29_Iprime;
    out.eq29.Q_pressure = Q29_pressure;
    out.eq29.Q_curvature = Q29_curvature;
    out.eq29.Q_poloidalField = Q29_poloidalField;
    out.eq29.Q = Q29;
    out.eq29.integrandJ = integrandJ;
    out.eq29.J = J;
    out.eq29.Jperiod = Jperiod;
    out.eq29.gradS2_metric = gradS2_metric;
    out.eq29.gradS2_shear = gradS2_shear;
    out.eq29.gradS2 = gradS2;

    out.eq31 = struct();
    out.eq31.toroidicity = curv_toroidicity;
    out.eq31.normalCurvature = curv_normal;
    out.eq31.magneticShear = curv_shear;
    out.eq31.bracket = curvBracket;
    out.eq31.drive = curvDrive;

    out.eq32 = struct();
    out.eq32.outerBp = Bp;
    out.eq32.A_inner = A_inner;
    out.eq32.C_drive = C_drive;
    out.eq32.C_dividedByBp = C_dividedByBp;
    out.eq32.description = ...
        'Bp*d/dl(A_inner*dFdl) + C_drive*F = 0';
    out.eq32.dividedDescription = ...
        'd/dl(A_inner*dFdl) + C_dividedByBp*F = 0';
end

%======================================================================
function ext = extend_ch3_terms(one, nside, opt)
    turns = (-nside):(nside-1);
    n = numel(one.l);
    nt = numel(turns);
    total = n * nt;

    ext = struct();
    ext.l = zeros(total, 1);
    ext.theta = zeros(total, 1);
    ext.turn = zeros(total, 1);

    names = {'R','Z','u','h0','Bp','Rc','Btor','B2','B4','dB2dl'};
    for i = 1:numel(names)
        ext.(names{i}) = zeros(total, 1);
    end

    eq29Names = fieldnames(one.eq29);
    eq31Names = fieldnames(one.eq31);
    eq32Names = {'outerBp','A_inner','C_drive','C_dividedByBp'};

    for i = 1:numel(eq29Names)
        value = one.eq29.(eq29Names{i});
        if isnumeric(value) && isvector(value) && numel(value) == n
            ext.eq29.(eq29Names{i}) = zeros(total, 1);
        else
            ext.eq29.(eq29Names{i}) = value;
        end
    end

    for i = 1:numel(eq31Names)
        value = one.eq31.(eq31Names{i});
        if isnumeric(value) && isvector(value) && numel(value) == n
            ext.eq31.(eq31Names{i}) = zeros(total, 1);
        else
            ext.eq31.(eq31Names{i}) = value;
        end
    end

    for i = 1:numel(eq32Names)
        ext.eq32.(eq32Names{i}) = zeros(total, 1);
    end

    cursor = 1;
    for k = turns
        idx = cursor:(cursor+n-1);

        ext.l(idx) = one.l + k * one.L;
        ext.theta(idx) = one.theta + 2*pi*k;
        ext.turn(idx) = k;

        ext.R(idx) = one.R;
        ext.Z(idx) = one.Z;
        ext.u(idx) = one.u - 2*pi*k;
        ext.h0(idx) = one.h0;
        ext.Bp(idx) = one.Bp;
        ext.Rc(idx) = one.Rc;
        ext.Btor(idx) = one.Btor;
        ext.B2(idx) = one.B2;
        ext.B4(idx) = one.B4;
        ext.dB2dl(idx) = one.dB2dl;

        % Periodic pieces repeat, but the ballooning integral continues
        % secularly from period to period.
        J = one.eq29.J + k * one.eq29.Jperiod;

        ext = copy_eq29_periodic(ext, one, idx);
        ext.eq29.J(idx) = J;

        gradS2_shear = (one.h0.^2 .* one.Bp.^2 .* one.I0^2 ./ one.x0^2) .* J.^2;
        gradS2 = one.eq29.gradS2_metric + gradS2_shear;

        ext.eq29.gradS2_shear(idx) = gradS2_shear;
        ext.eq29.gradS2(idx) = gradS2;

        curv_shear = -one.I0^2 .* one.h0 .* one.Bp ./ ...
            (one.x0 .* one.B4) .* one.dB2dl .* J;
        curvBracket = one.eq31.toroidicity + one.eq31.normalCurvature + curv_shear;
        curvDrive = opt.DriveFactor .* one.pprimeEquation ./ ...
            (one.x0 .* one.h0) .* curvBracket;

        ext.eq31.toroidicity(idx) = one.eq31.toroidicity;
        ext.eq31.normalCurvature(idx) = one.eq31.normalCurvature;
        ext.eq31.magneticShear(idx) = curv_shear;
        ext.eq31.bracket(idx) = curvBracket;
        ext.eq31.drive(idx) = curvDrive;

        A_inner = gradS2 ./ one.B2 .* one.Bp;
        ext.eq32.outerBp(idx) = one.Bp;
        ext.eq32.A_inner(idx) = A_inner;
        ext.eq32.C_drive(idx) = curvDrive;
        ext.eq32.C_dividedByBp(idx) = curvDrive ./ one.Bp;

        cursor = cursor + n;
    end

    ext.L = one.L;
    ext.x0 = one.x0;
    ext.I0 = one.I0;
    ext.Iprime0 = one.Iprime0;
    ext.pprimeEquation = one.pprimeEquation;
    ext.direction = one.direction;
    ext.NPeriodsEachSide = nside;
    ext.eq29.Jperiod = one.eq29.Jperiod;
    ext.eq32.description = one.eq32.description;
    ext.eq32.dividedDescription = one.eq32.dividedDescription;

    if ~isempty(opt.ModeF)
        ext.modeTerms = evaluate_mode_terms(ext, opt.ModeF(:));
    else
        ext.modeTerms = [];
    end
end

%======================================================================
function ext = copy_eq29_periodic(ext, one, idx)
    names = fieldnames(one.eq29);
    for i = 1:numel(names)
        name = names{i};
        if any(strcmp(name, {'J','gradS2_shear','gradS2','Jperiod'}))
            continue;
        end

        value = one.eq29.(name);
        if isnumeric(value) && isvector(value) && numel(value) == numel(one.l)
            ext.eq29.(name)(idx) = value;
        end
    end
end

%======================================================================
function mode = evaluate_mode_terms(path, F)
    if numel(F) ~= numel(path.theta)
        error('compute_bishop_ch3_terms:ModeFSizeMismatch', ...
            'ModeF must have the same length as the selected theta grid.');
    end

    l = path.l(:);
    F = F(:);
    A = path.eq32.A_inner(:);
    Bp = path.eq32.outerBp(:);
    C = path.eq32.C_drive(:);

    dFdl = gradient(F, l);
    flux = A .* dFdl;
    selfAdjointBending = gradient(flux, l);
    fieldLineBending = Bp .* selfAdjointBending;
    pressureCurvature = C .* F;

    mode = struct();
    mode.F = F;
    mode.dFdl = dFdl;
    mode.flux = flux;
    mode.selfAdjointBending = selfAdjointBending;
    mode.fieldLineBending = fieldLineBending;
    mode.pressureCurvature = pressureCurvature;
    mode.residual = fieldLineBending + pressureCurvature;
    mode.dividedResidual = selfAdjointBending + path.eq32.C_dividedByBp(:).*F;
end

%======================================================================
function [J, Jperiod] = periodic_cumtrapz(l, f, L)
    l = l(:);
    f = f(:);

    lEnd = [l; L];
    fEnd = [f; f(1)];
    dl = diff(lEnd);

    JEnd = [0; cumsum(0.5 .* (fEnd(1:end-1) + fEnd(2:end)) .* dl)];
    J = JEnd(1:end-1);
    Jperiod = JEnd(end);
end

%======================================================================
function dydl = periodic_derivative(l, y, L)
    l = l(:);
    y = y(:);
    n = numel(l);

    if n < 3
        error('compute_bishop_ch3_terms:TooFewPoints', ...
            'Need at least three points for derivatives.');
    end

    dlForward = [diff(l); L - l(end) + l(1)];
    dydl = zeros(n, 1);

    for i = 1:n
        im = i - 1; if im < 1, im = n; end
        ip = i + 1; if ip > n, ip = 1; end

        ym = y(im);
        yp = y(ip);
        if ip == 1
            % y is periodic, but l(ip) is one period ahead of l(i).
            yp = y(1);
        end
        if im == n
            ym = y(n);
        end

        ds = dlForward(im) + dlForward(i);
        dydl(i) = (yp - ym) / ds;
    end
end

%======================================================================
function plot_ch3_terms(ch3)
    figure('Color', 'w', 'Name', 'Bishop Chapter 3 terms');
    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(ch3.theta, ch3.Bp, 'LineWidth', 1.3);
    grid on;
    xlabel('theta');
    ylabel('Bp');
    title('Poloidal field');

    nexttile;
    plot(ch3.theta, ch3.eq29.J, 'LineWidth', 1.3);
    grid on;
    xlabel('theta');
    ylabel('J');
    title('Eq. 29 integral');

    nexttile;
    plot(ch3.theta, ch3.eq29.gradS2, 'LineWidth', 1.3);
    hold on;
    plot(ch3.theta, ch3.eq29.gradS2_metric, '--', 'LineWidth', 1.0);
    plot(ch3.theta, ch3.eq29.gradS2_shear, ':', 'LineWidth', 1.0);
    grid on;
    xlabel('theta');
    ylabel('|grad S|^2');
    title('Eq. 29');
    legend('total', 'metric', 'shear', 'Location', 'best');

    nexttile;
    plot(ch3.theta, ch3.eq31.drive, 'LineWidth', 1.3);
    grid on;
    xlabel('theta');
    ylabel('C');
    title('Eq. 31 drive');

    nexttile;
    plot(ch3.theta, ch3.eq32.A_inner, 'LineWidth', 1.3);
    grid on;
    xlabel('theta');
    ylabel('A inner');
    title('Field-line bending coefficient');

    nexttile;
    plot(ch3.theta, ch3.eq32.C_dividedByBp, 'LineWidth', 1.3);
    grid on;
    xlabel('theta');
    ylabel('C/Bp');
    title('Divided drive coefficient');

    if isfield(ch3, 'extended') && ~isempty(ch3.extended)
        figure('Color', 'w', 'Name', 'Bishop Chapter 3 extended theta');
        tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

        ex = ch3.extended;

        nexttile;
        plot(ex.theta, ex.eq29.J, 'LineWidth', 1.2);
        grid on;
        xlabel('theta');
        ylabel('J');
        title('Secular ballooning integral');

        nexttile;
        plot(ex.theta, ex.eq29.gradS2, 'LineWidth', 1.2);
        grid on;
        xlabel('theta');
        ylabel('|grad S|^2');
        title('Extended Eq. 29');

        nexttile;
        plot(ex.theta, ex.eq31.drive, 'LineWidth', 1.2);
        grid on;
        xlabel('theta');
        ylabel('C');
        title('Extended Eq. 31 drive');

        nexttile;
        plot(ex.theta, ex.eq32.A_inner, 'LineWidth', 1.2);
        grid on;
        xlabel('theta');
        ylabel('A inner');
        title('Extended bending coefficient');
    end
end
