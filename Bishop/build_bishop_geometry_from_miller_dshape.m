function geom = build_bishop_geometry_from_miller_dshape(varargin)
%BUILD_BISHOP_GEOMETRY_FROM_MILLER_DSHAPE Convert Miller D-shape data to Bishop geom.
%
%   geom = BUILD_BISHOP_GEOMETRY_FROM_MILLER_DSHAPE()
%   geom = BUILD_BISHOP_GEOMETRY_FROM_MILLER_DSHAPE('Name', value, ...)
%
% This adapter uses the D-shape Miller example in ../Miller and returns a
% geometry structure compatible with COMPUTE_BISHOP_CH3_TERMS and the
% s-alpha scan utilities. The Miller curve is re-ordered to Bishop's
% clockwise l convention and resampled uniformly in arclength.
%
% The nominal Miller point is normalized with Miller's alpha definition,
%
%   alpha = -2*q^2*R0/B0^2 * dp/dr
%         = -2*q^2*R0*(dpsi/dr)/B0^2 * dp/dpsi,
%
% and by enforcing
%
%   Jperiod = -2*pi*s_hat
%
% for the baseline p.alpha and p.s_hat. This gives the s-alpha scan a
% canonical Miller-style shear convention while preserving Bishop's Eq. (29)
% construction.

    ip = inputParser;
    addParameter(ip, 'MillerPath', fullfile(fileparts(pwd), 'Miller'), @(x)ischar(x) || isstring(x));
    addParameter(ip, 'Params', [], @(x)isstruct(x) || isempty(x));
    addParameter(ip, 'SourceNTheta', 801, @(x)isnumeric(x) && isscalar(x) && x >= 16);
    addParameter(ip, 'NTheta', 256, @(x)isnumeric(x) && isscalar(x) && x >= 16);
    addParameter(ip, 'JShearFactor', -2*pi, @(x)isnumeric(x) && isscalar(x));
    parse(ip, varargin{:});
    opt = ip.Results;

    millerPath = char(opt.MillerPath);
    if exist(millerPath, 'dir') ~= 7
        error('build_bishop_geometry_from_miller_dshape:MissingMillerPath', ...
            'MillerPath does not exist: %s', millerPath);
    end
    addpath(millerPath);

    if isempty(opt.Params)
        p = DshapeMillerParams();
    else
        p = opt.Params;
    end

    [BpSource, out] = compute_miller_bpol_corrected(p, round(opt.SourceNTheta));
    [R, Z, Bp, ~, ~] = resample_miller_curve_clockwise( ...
        out.R(:), out.Z(:), BpSource(:), round(opt.NTheta));

    tn = compute_tangent_normal(R, Z, ...
        'Closed', true, ...
        'Method', 'central', ...
        'NormalDirection', 'bishop', ...
        'Direction', 'clockwise', ...
        'RemoveDuplicateLastPoint', true);
    curv = compute_curvature(tn.R, tn.Z, ...
        'Closed', true, ...
        'Method', 'central', ...
        'NormalDirection', 'bishop', ...
        'Direction', 'clockwise', ...
        'RemoveDuplicateLastPoint', false);

    R = tn.R(:);
    Z = tn.Z(:);
    l = tn.l(:);
    L = tn.L;
    Bp = max(abs(Bp(:)), 1e-12);
    if numel(Bp) ~= numel(R)
        error('build_bishop_geometry_from_miller_dshape:BpSizeMismatch', ...
            'Resampled Bp and geometry grids have different lengths.');
    end

    R0 = p.A .* p.r;
    Z0 = 0;
    I0 = out.F;
    q0 = p.q;
    Btor = I0 ./ R;
    B2 = Bp.^2 + Btor.^2;
    gradpsi = R .* Bp;

    alphaFactor = -2 .* q0.^2 .* R0 .* out.dpdr ./ (p.B0.^2);
    mu0Pprime0 = p.alpha ./ alphaFactor;

    geom = struct();
    geom.R = R;
    geom.Z = Z;
    geom.l = l;
    geom.L = L;
    geom.u = tn.u(:);
    geom.h0 = R ./ R(1);
    geom.Bp = Bp(:);
    geom.gradpsi = gradpsi(:);
    geom.R0 = R0;
    geom.Z0 = Z0;
    geom.x0 = R(1);
    geom.I0 = I0;
    geom.Iprime0 = 0;
    geom.pprime0 = mu0Pprime0 ./ (4*pi*1e-7);
    geom.mu0_pprime0 = mu0Pprime0;
    geom.q0 = q0;
    geom.shear0 = p.s_hat;
    geom.psi_axis = 0;
    geom.psi_target = 0.5;
    geom.psi_relative = 0.5;
    geom.direction = 'clockwise';
    geom.source = 'Miller D-shape local equilibrium';

    copyNames = {'tR','tZ','nR','nZ','dRdl','dZdl','dl','closed'};
    for i = 1:numel(copyNames)
        name = copyNames{i};
        if isfield(tn, name)
            geom.(name) = tn.(name);
        end
    end

    curvatureNames = {'kappa_signed','kappa','Rc_signed','Rc','dudl','dtRdl','dtZdl'};
    for i = 1:numel(curvatureNames)
        name = curvatureNames{i};
        geom.(name) = curv.(name);
    end

    geom.Btor = Btor(:);
    geom.B2 = B2(:);
    geom.alphaFactor0 = alphaFactor;
    geom.JShearFactor0 = opt.JShearFactor;
    geom.shape = estimate_miller_shape(R, Z, R0, Z0);
    geom.miller = struct();
    geom.miller.params = p;
    geom.miller.output = out;
    geom.miller.sourceNTheta = round(opt.SourceNTheta);
    geom.miller.NTheta = round(opt.NTheta);
    geom.miller.dpdr_flux = out.dpdr;
    geom.miller.alphaFactor = alphaFactor;
    geom.miller.q_check = out.q_check;

    geom.Iprime0 = baseline_iprime_for_target_jperiod(geom, opt.JShearFactor .* p.s_hat);
end

%======================================================================
function [Bp, out] = compute_miller_bpol_corrected(p, ntheta)
    [Ru, Zu, Rr, Zr, jac, u, R, Z, p] = DshapeDerivatives(p, ntheta);

    R0 = p.A .* p.r;
    if isfield(p, 'F')
        F = p.F;
    else
        F = R0 .* p.B0;
    end

    hu = hypot(Ru, Zu);
    if isfield(p, 'dpdr')
        dpdr = p.dpdr;
    else
        qIntegral = trapz(u, abs(jac) ./ (R .* hu));
        dpdr = F .* p.r .* qIntegral ./ (2*pi*p.q);
    end

    Bp = abs(dpdr) .* hu ./ (abs(jac) .* R);

    out = struct();
    out.u = u;
    out.R = R;
    out.Z = Z;
    out.Ru = Ru;
    out.Zu = Zu;
    out.Rr = Rr;
    out.Zr = Zr;
    out.jac = jac;
    out.hu = hu;
    out.F = F;
    out.dpdr = dpdr;
    out.dpsidr = dpdr;
    out.q_integral = trapz(u, abs(jac) ./ (R .* hu));
    out.q_check = F .* p.r .* out.q_integral ./ (2*pi*dpdr);
    out.params = p;
end

%======================================================================
function [R, Z, Bp, l, L] = resample_miller_curve_clockwise(R0, Z0, Bp0, ntheta)
    R0 = R0(:);
    Z0 = Z0(:);
    Bp0 = Bp0(:);
    if numel(R0) ~= numel(Z0) || numel(R0) ~= numel(Bp0)
        error('build_bishop_geometry_from_miller_dshape:BadSourceSize', ...
            'Source R, Z, and Bp arrays must have the same length.');
    end

    if hypot(R0(end) - R0(1), Z0(end) - Z0(1)) < 1e-10
        R0 = R0(1:end-1);
        Z0 = Z0(1:end-1);
        Bp0 = Bp0(1:end-1);
    end

    % Miller's u increases counter-clockwise from the outboard midplane.
    % Bishop's report draws l clockwise, so reverse the non-initial points.
    idx = [1; (numel(R0):-1:2).'];

    Rsrc = R0(idx);
    Zsrc = Z0(idx);
    Bpsrc = Bp0(idx);

    dseg = hypot(diff([Rsrc; Rsrc(1)]), diff([Zsrc; Zsrc(1)]));
    L = sum(dseg);
    lNodes = [0; cumsum(dseg(1:end-1))];
    lExt = [lNodes; L];

    l = linspace(0, L, ntheta + 1).';
    l(end) = [];

    R = interp1(lExt, [Rsrc; Rsrc(1)], l, 'pchip');
    Z = interp1(lExt, [Zsrc; Zsrc(1)], l, 'pchip');
    Bp = interp1(lExt, [Bpsrc; Bpsrc(1)], l, 'pchip');
end

%======================================================================
function Iprime0 = baseline_iprime_for_target_jperiod(geom, targetJperiod)
    opt = {'L0Index', 1, 'NPeriodsEachSide', 0, 'Plot', false};

    J0 = evaluate_jperiod_local(geom, 0, 0, opt);
    JI = evaluate_jperiod_local(geom, 1, 0, opt) - J0;
    Jp = evaluate_jperiod_local(geom, 0, 1, opt) - J0;

    if abs(JI) < 1e-12
        error('build_bishop_geometry_from_miller_dshape:WeakIprimeControl', ...
            'Bishop Eq. (29) Jperiod is insensitive to Iprime.');
    end

    Iprime0 = (targetJperiod - J0 - Jp .* geom.mu0_pprime0) ./ JI;
end

%======================================================================
function Jperiod = evaluate_jperiod_local(geom, Iprime0, pprimeEquation, opt)
    ch3 = compute_bishop_ch3_terms(geom, ...
        opt{:}, ...
        'I0', geom.I0, ...
        'Iprime0', Iprime0, ...
        'PprimeEquation', pprimeEquation);
    Jperiod = ch3.eq29.Jperiod;
end

%======================================================================
function shape = estimate_miller_shape(R, Z, R0, Z0)
    R = R(:);
    Z = Z(:);
    [Rmax, imax] = max(R);
    Rmin = min(R);
    [Zmax, itop] = max(Z);
    [Zmin, ibot] = min(Z);
    a = 0.5 .* (Rmax - Rmin);

    shape = struct();
    shape.Raxis = R0;
    shape.Zaxis = Z0;
    shape.Rmax = Rmax;
    shape.Rmin = Rmin;
    shape.Zmax = Zmax;
    shape.Zmin = Zmin;
    shape.minorRadius = a;
    shape.elongation = 0.5 .* (Zmax - Zmin) ./ max(a, eps);
    shape.epsilon = a ./ R0;
    shape.deltaTop = (R0 - R(itop)) ./ max(a, eps);
    shape.deltaBottom = (R0 - R(ibot)) ./ max(a, eps);
    shape.deltaMean = 0.5 .* (shape.deltaTop + shape.deltaBottom);
    shape.outboardIndex = imax;
    shape.topIndex = itop;
    shape.bottomIndex = ibot;
end
