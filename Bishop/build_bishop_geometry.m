function geom = build_bishop_geometry(eq, psin_target, varargin)
%BUILD_BISHOP_GEOMETRY Prepare local geometry and flux-function data
% for a chosen flux surface.

    p = inputParser;
    addParameter(p, 'UseUniformTheta', false, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'UseUniformArclength', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'NTheta', 256, @(x)isnumeric(x) && isscalar(x) && x >= 8);
    addParameter(p, 'ThetaCenter', 'axis', @(x)ischar(x) || isstring(x));
    addParameter(p, 'InterpMethod', 'linear', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Direction', 'clockwise', @(x)ischar(x) || isstring(x));
    addParameter(p, 'NormalDirection', 'bishop', @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    opt = p.Results;

    mu0 = 4*pi*1e-7;

    %-----------------------------------
    % 1) Extract flux surface
    %-----------------------------------
    direction = lower(string(opt.Direction));
    if direction == "clockwise"
        clockwise = true;
    elseif direction == "counterclockwise"
        clockwise = false;
    else
        clockwise = [];
    end

    surf = extract_flux_surface(eq, psin_target, ...
        'Clockwise', clockwise, ...
        'UniformTheta', opt.UseUniformTheta, ...
        'NTheta', opt.NTheta, ...
        'ThetaCenter', opt.ThetaCenter);

    if opt.UseUniformArclength
        surf_l = resample_closed_curve_arclength(surf.R, surf.Z, opt.NTheta, ...
            'Direction', opt.Direction);
        surf.R = surf_l.R;
        surf.Z = surf_l.Z;
    end

    %-----------------------------------
    % 2) Arc length, tangent, Bishop normal, and u
    %-----------------------------------
    tn = compute_tangent_normal(surf.R, surf.Z, ...
        'Closed', true, ...
        'Method', 'central', ...
        'NormalDirection', opt.NormalDirection, ...
        'Direction', opt.Direction);

    R = tn.R(:);
    Z = tn.Z(:);
    l = tn.l(:);
    u = unwrap(atan2(tn.tZ, tn.tR));

    %-----------------------------------
    % 3) Curvature
    %-----------------------------------
    curv = compute_curvature(R, Z, ...
        'Closed', true, ...
        'Method', 'central', ...
        'NormalDirection', opt.NormalDirection, ...
        'Direction', opt.Direction);

    kappa_signed = curv.kappa_signed;
    kappa = curv.kappa;
    Rc_signed = curv.Rc_signed;
    Rc = curv.Rc;

    %-----------------------------------
    % 4) h0 = R / R0
    %-----------------------------------
    R0 = R(1);
    h0 = R / R0;

    %-----------------------------------
    % 5) psi interpolation and grad(psi)
    %-----------------------------------
    [Rgrid, Zgrid, psiRZ] = local_grids(eq);

    psi_fun = griddedInterpolant({Rgrid, Zgrid}, psiRZ, ...
        char(opt.InterpMethod), 'none');

    [dpsi_dZ_grid, dpsi_dR_grid] = gradient(psiRZ, Zgrid, Rgrid);

    dpsiR_fun = griddedInterpolant({Rgrid, Zgrid}, dpsi_dR_grid, ...
        char(opt.InterpMethod), 'none');
    dpsiZ_fun = griddedInterpolant({Rgrid, Zgrid}, dpsi_dZ_grid, ...
        char(opt.InterpMethod), 'none');

    psi_vals = psi_fun(R, Z);
    dpsiR = dpsiR_fun(R, Z);
    dpsiZ = dpsiZ_fun(R, Z);

    gradpsi = hypot(dpsiR, dpsiZ);
    Bp = gradpsi ./ R;
    shape = estimate_shape(R, Z, eq.rmaxis, eq.zmaxis);

    %-----------------------------------
    % 6) Flux-function profiles.
    %-----------------------------------
    % GEQDSK profiles are tabulated versus the usual outward-increasing
    % poloidal flux. Bishop's equations use this flux derivative directly;
    % the usual positive s-alpha pressure-gradient coordinate is introduced
    % later by alpha = -const*(mu0*dp/dpsi).
    psin_grid = linspace(0, 1, eq.nw);
    interpMethod = char(opt.InterpMethod);

    I0      = interp1(psin_grid, eq.fpol,   psin_target, interpMethod, 'extrap');
    ffprim0 = interp1(psin_grid, eq.ffprim, psin_target, interpMethod, 'extrap');
    pprime0 = interp1(psin_grid, eq.pprime, psin_target, interpMethod, 'extrap');
    q0      = interp1(psin_grid, eq.qpsi,   psin_target, interpMethod, 'extrap');

    if numel(eq.qpsi) >= 3
        dq_dpsin_grid = gradient(eq.qpsi(:), psin_grid(:));
        dq_dpsin0 = interp1(psin_grid, dq_dpsin_grid, psin_target, interpMethod, 'extrap');
    else
        dq_dpsin0 = NaN;
    end

    mu0_pprime_geqdsk0 = mu0 * pprime0;
    dpsi_total = eq.sibry - eq.simag;
    if abs(dpsi_total) > 0
        dq_dpsi0 = dq_dpsin0 ./ dpsi_total;
    else
        dq_dpsi0 = NaN;
    end

    psi_target = eq.simag + psin_target * dpsi_total;
    psi_relative = psi_target - eq.simag;
    if abs(q0) > 0
        shear0 = 2 * psi_relative .* dq_dpsi0 ./ q0;
    else
        shear0 = NaN;
    end

    if abs(I0) < 1e-12
        warning('build_bishop_geometry:InearZero', ...
            'I is near zero on this surface; setting Iprime = NaN.');
        Iprime_geqdsk0 = NaN;
    else
        Iprime_geqdsk0 = ffprim0 / I0;
    end

    I          = I0                  * ones(size(l));
    Iprime     = Iprime_geqdsk0      * ones(size(l));
    pprime     = pprime0             * ones(size(l));
    mu0_pprime = mu0_pprime_geqdsk0  * ones(size(l));
    ffprim     = ffprim0             * ones(size(l));

    %-----------------------------------
    % 8) Output
    %-----------------------------------
    geom = struct();

    geom.R = R;
    geom.Z = Z;
    geom.l = l;
    geom.L = tn.L;
    geom.dl = tn.dl;

    geom.tR = tn.tR(:);
    geom.tZ = tn.tZ(:);
    geom.nR = tn.nR(:);
    geom.nZ = tn.nZ(:);

    geom.u = u;

    geom.kappa_signed = kappa_signed(:);
    geom.kappa = kappa(:);
    geom.Rc_signed = Rc_signed(:);
    geom.Rc = Rc(:);

    geom.h0 = h0(:);

    geom.psi = psi_vals(:);
    geom.psin = psin_target * ones(size(l));
    geom.psi_axis = eq.simag;
    geom.psi_bdry = eq.sibry;
    geom.psi_target = psi_target;
    geom.psi_relative = psi_relative;

    geom.dpsi_dR = dpsiR(:);
    geom.dpsi_dZ = dpsiZ(:);
    geom.gradpsi = gradpsi(:);
    geom.Bp = Bp(:);

    geom.R0 = R0;
    geom.direction = tn.direction;
    geom.normalDirection = tn.normalDirection;
    geom.theta_ballooning = 2*pi * l(:) / tn.L;
    geom.shape = shape;

    % flux-function scalars
    geom.I0 = I0;
    geom.Iprime0 = Iprime_geqdsk0;
    geom.pprime0 = pprime0;
    geom.mu0_pprime0 = mu0_pprime_geqdsk0;
    geom.ffprim0 = ffprim0;
    geom.Iprime_geqdsk0 = Iprime_geqdsk0;
    geom.pprime_geqdsk0 = pprime0;
    geom.mu0_pprime_geqdsk0 = mu0_pprime_geqdsk0;
    geom.ffprim_geqdsk0 = ffprim0;
    geom.q0 = q0;
    geom.dq_dpsin0 = dq_dpsin0;
    geom.dq_dpsi0 = dq_dpsi0;
    geom.shear0 = shear0;

    % same data mapped onto l
    geom.I = I(:);
    geom.Iprime = Iprime(:);
    geom.pprime = pprime(:);
    geom.mu0_pprime = mu0_pprime(:);
    geom.ffprim = ffprim(:);
end

% -------------------------------------------------------------------------
function shape = estimate_shape(R, Z, Raxis, Zaxis)
    R = R(:);
    Z = Z(:);

    Rmax = max(R);
    Rmin = min(R);
    Zmax = max(Z);
    Zmin = min(Z);

    a = 0.5 * (Rmax - Rmin);
    b = 0.5 * (Zmax - Zmin);

    [~, itop] = max(Z);
    [~, ibot] = min(Z);

    shape = struct();
    shape.Raxis = Raxis;
    shape.Zaxis = Zaxis;
    shape.Rmax = Rmax;
    shape.Rmin = Rmin;
    shape.Zmax = Zmax;
    shape.Zmin = Zmin;
    shape.minorRadius = a;
    shape.elongation = b ./ a;
    shape.epsilon = a ./ Raxis;
    shape.deltaTop = (Raxis - R(itop)) ./ a;
    shape.deltaBottom = (Raxis - R(ibot)) ./ a;
    shape.deltaMean = 0.5 .* (shape.deltaTop + shape.deltaBottom);
end

% -------------------------------------------------------------------------
function [Rgrid, Zgrid, psiRZ] = local_grids(eq)
    Rgrid = eq.rgrid(:);
    Zgrid = eq.zgrid(:);

    psi = eq.psirz;
    if isequal(size(psi), [numel(Rgrid), numel(Zgrid)])
        psiRZ = psi;
    elseif isequal(size(psi), [numel(Zgrid), numel(Rgrid)])
        psiRZ = psi.';
    else
        error('build_bishop_geometry:BadPsiGrid', ...
            'eq.psirz size does not match R/Z grid.');
    end
end
