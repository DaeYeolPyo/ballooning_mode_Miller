function cont = eq_straight_fieldline_theta(eq, psin_target, varargin)
%EQ_STRAIGHT_FIELDLINE_THETA Compute flux-surface theta_coord from GEQDSK equilibrium.
%
%   cont = eq_straight_fieldline_theta(eq, psin_target)
%   cont = eq_straight_fieldline_theta(eq, psin_target, 'Name', value, ...)
%
% This is a MATLAB translation of the main logic in eq_coordinate.f90:
%   1) construct normalized psi map
%   2) find points on a target flux surface by radial bisection from axis
%      to the LCFS/boundary points
%   3) compute grad(psi), Bpol, and a coordinate-dependent Jacobian
%   4) integrate R/(Jacobian*|grad(psi)|) along the contour
%   5) normalize theta_coord to [0, 2*pi)
%
% Required eq fields, using common GEQDSK names:
%   rdim, zdim, rleft, zmid, rmaxis, zmaxis, simag, sibry, psirz
%   rbbbs, zbbbs
%
% Optional eq fields:
%   rgrid, zgrid, nw, nh
%
% Name-value options:
%   'JacMode'      : 1 equal-arc-length, 2 Hamada, 3 Boozer, 4 PEST. Default 4.
%   'NTheta'       : number of output points. Default length(eq.rbbbs).
%   'SortContour'  : true/false, sort by geometric angle. Default true.
%   'Tol'          : bisection tolerance. Default 1e-8.
%   'MaxIter'      : bisection max iterations. Default 100.
%   'InterpMethod' : 'spline', 'makima', or 'linear'. Default 'spline'.
%
% Output fields:
%   cont.R, cont.Z, cont.theta_geo, cont.theta_coord, cont.psin, cont.psi,
%   cont.JacMode, cont.integrand

    p = inputParser;
    addParameter(p, 'JacMode', 4, @(x)isnumeric(x) && isscalar(x) && any(x == 1:4));
    addParameter(p, 'NTheta', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 8));
    addParameter(p, 'SortContour', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Tol', 1e-8, @(x)isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'MaxIter', 100, @(x)isnumeric(x) && isscalar(x) && x >= 10);
    addParameter(p, 'InterpMethod', 'spline', @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    opt = p.Results;
    opt.InterpMethod = char(lower(string(opt.InterpMethod)));
    if strcmp(opt.InterpMethod, 'pchip')
        opt.InterpMethod = 'spline';   % griddedInterpolant does not support pchip in 2D
    end

    if psin_target < 0 || psin_target > 1
        error('psin_target must be in [0, 1].');
    end

    [Rgrid, Zgrid, psiRZ] = local_grids(eq);
    psiN = (psiRZ - eq.simag) ./ (eq.sibry - eq.simag);
    psi_target = eq.simag + psin_target * (eq.sibry - eq.simag);

    psiF = griddedInterpolant({Rgrid, Zgrid}, psiN, opt.InterpMethod, 'none');

    [gradPsi, Bpol] = local_grad_bpol(eq, Rgrid, Zgrid, psiRZ);
    Jacobian = local_jacobian(opt.JacMode, Rgrid, gradPsi, Bpol);

    gradF = griddedInterpolant({Rgrid, Zgrid}, gradPsi, opt.InterpMethod, 'none');
    jacF  = griddedInterpolant({Rgrid, Zgrid}, Jacobian, opt.InterpMethod, 'none');

    Rsep = eq.rbbbs(:);
    Zsep = eq.zbbbs(:);

    if isempty(opt.NTheta)
        nout = numel(Rsep);
    else
        nout = round(opt.NTheta);
    end

    % Use LCFS points as radial directions, then optionally resample directions
    % to a uniform geometric theta grid.
    theta_sep = mod(atan2(Zsep - eq.zmaxis, Rsep - eq.rmaxis), 2*pi);
    [theta_sep, iu] = unique(theta_sep, 'stable');
    Rsep = Rsep(iu);
    Zsep = Zsep(iu);

    if opt.SortContour
        [theta_sep, ord] = sort(theta_sep);
        Rsep = Rsep(ord);
        Zsep = Zsep(ord);
    end

    theta_dir = linspace(0, 2*pi, nout + 1).';
    theta_dir(end) = [];

    Rsep_ext = [Rsep; Rsep; Rsep];
    Zsep_ext = [Zsep; Zsep; Zsep];
    theta_ext = [theta_sep - 2*pi; theta_sep; theta_sep + 2*pi];

    % MATLAB interp1 requires unique sample points.
    % LCFS data can contain repeated or nearly repeated poloidal angles,
    % especially around the 0/2pi branch cut. Clean them before interp1.
    [theta_ext, Rsep_ext, Zsep_ext] = unique_theta_samples(theta_ext, Rsep_ext, Zsep_ext);

    Rdir = interp1(theta_ext, Rsep_ext, theta_dir, 'pchip');
    Zdir = interp1(theta_ext, Zsep_ext, theta_dir, 'pchip');

    R = zeros(nout, 1);
    Z = zeros(nout, 1);

    if abs(psin_target - 1) < opt.Tol
        R = Rdir;
        Z = Zdir;
    else
        for i = 1:nout
            [R(i), Z(i)] = radial_bisection_to_flux( ...
                psiF, psin_target, eq.rmaxis, eq.zmaxis, Rdir(i), Zdir(i), ...
                opt.Tol, opt.MaxIter);
        end
    end

    theta_geo = mod(atan2(Z - eq.zmaxis, R - eq.rmaxis), 2*pi);
    if opt.SortContour
        [theta_geo, ord] = sort(theta_geo);
        R = R(ord);
        Z = Z(ord);
    end

    [theta_coord, integrand] = evaluate_poloidal_angle(R, Z, gradF, jacF);

    cont = struct();
    cont.psin = psin_target;
    cont.psi = psi_target;
    cont.R = R(:);
    cont.Z = Z(:);
    cont.theta_geo = theta_geo(:);
    cont.theta_coord = theta_coord(:);
    cont.JacMode = opt.JacMode;
    cont.integrand = integrand(:);
end

% -------------------------------------------------------------------------
function [theta_u, R_u, Z_u] = unique_theta_samples(theta, R, Z)
    theta = theta(:);
    R = R(:);
    Z = Z(:);

    good = isfinite(theta) & isfinite(R) & isfinite(Z);
    theta = theta(good);
    R = R(good);
    Z = Z(good);

    [theta, ord] = sort(theta);
    R = R(ord);
    Z = Z(ord);

    if isempty(theta)
        theta_u = theta;
        R_u = R;
        Z_u = Z;
        return;
    end

    tol = 1e-12 * max(1, max(abs(theta)));
    group = cumsum([1; abs(diff(theta)) > tol]);
    ng = group(end);

    theta_u = zeros(ng, 1);
    R_u = zeros(ng, 1);
    Z_u = zeros(ng, 1);

    for k = 1:ng
        m = (group == k);
        theta_u(k) = mean(theta(m));
        R_u(k) = mean(R(m));
        Z_u(k) = mean(Z(m));
    end
end

% -------------------------------------------------------------------------
function [Rgrid, Zgrid, psiRZ] = local_grids(eq)
    if isfield(eq, 'rgrid')
        Rgrid = eq.rgrid(:);
    else
        nw = local_get_nw(eq);
        Rgrid = eq.rleft + (0:nw-1).' * eq.rdim / (nw - 1);
    end

    if isfield(eq, 'zgrid')
        Zgrid = eq.zgrid(:);
    else
        nh = local_get_nh(eq);
        Zgrid = eq.zmid - eq.zdim/2 + (0:nh-1).' * eq.zdim / (nh - 1);
    end

    psi = eq.psirz;
    if isequal(size(psi), [numel(Rgrid), numel(Zgrid)])
        psiRZ = psi;
    elseif isequal(size(psi), [numel(Zgrid), numel(Rgrid)])
        psiRZ = psi.';
    else
        error('eq.psirz size does not match R/Z grid.');
    end
end

function nw = local_get_nw(eq)
    if isfield(eq, 'nw'), nw = eq.nw; else, nw = size(eq.psirz, 1); end
end

function nh = local_get_nh(eq)
    if isfield(eq, 'nh'), nh = eq.nh; else, nh = size(eq.psirz, 2); end
end

% -------------------------------------------------------------------------
function [gradPsi, Bpol] = local_grad_bpol(~, Rgrid, Zgrid, psiRZ)
    % psiRZ is nR x nZ. MATLAB gradient takes the column coordinate first
    % and the row coordinate second, so columns use Z and rows use R here.
    [psi_Z, psi_R] = gradient(psiRZ, Zgrid, Rgrid);

    RR = repmat(Rgrid(:), 1, numel(Zgrid));
    BR = -psi_Z ./ RR;
    BZ =  psi_R ./ RR;

    gradPsi = hypot(psi_R, psi_Z);
    Bpol = hypot(BR, BZ);
end

% -------------------------------------------------------------------------
function J = local_jacobian(jacMode, Rgrid, gradPsi, Bpol)
    RR = repmat(Rgrid(:), 1, size(gradPsi, 2));

    switch jacMode
        case 1
            % equal-arc-length mode in the original Fortran
            J = RR ./ gradPsi;
        case 2
            % Hamada
            J = ones(size(gradPsi));
        case 3
            % Boozer
            J = 1 ./ (Bpol.^2);
        case 4
            % PEST
            J = RR.^2;
    end

    J(~isfinite(J)) = NaN;
end

% -------------------------------------------------------------------------
function [Rroot, Zroot] = radial_bisection_to_flux(psiF, target, Raxis, Zaxis, Rsep, Zsep, tol, maxIter)
    aR = Raxis; aZ = Zaxis;
    bR = Rsep;  bZ = Zsep;

    fa = psiF(aR, aZ) - target;
    fb = psiF(bR, bZ) - target;

    if ~isfinite(fa) || ~isfinite(fb)
        error('Interpolation returned NaN at axis or boundary.');
    end

    if fa * fb > 0
        error('Bisection failed: target flux is not bracketed along one radial ray.');
    end

    cR = NaN; cZ = NaN;
    for it = 1:maxIter
        cR = 0.5 * (aR + bR);
        cZ = 0.5 * (aZ + bZ);
        fc = psiF(cR, cZ) - target;

        if abs(fc) < tol || hypot(bR-aR, bZ-aZ) < tol
            break;
        end

        if fa * fc < 0
            bR = cR; bZ = cZ; fb = fc; %#ok<NASGU>
        else
            aR = cR; aZ = cZ; fa = fc;
        end
    end

    Rroot = cR;
    Zroot = cZ;
end

% -------------------------------------------------------------------------
function [theta_coord, integrand] = evaluate_poloidal_angle(R, Z, gradF, jacF)
    R = R(:); Z = Z(:);
    n = numel(R);

    Rc = [R; R(1)];
    Zc = [Z; Z(1)];

    gradP = gradF(Rc, Zc);
    Jac   = jacF(Rc, Zc);

    f = Rc ./ (Jac .* gradP);
    dlp = hypot(diff(Rc), diff(Zc));

    seg = dlp .* 0.5 .* (f(1:end-1) + f(2:end));
    s = [0; cumsum(seg(1:end-1))];
    total = sum(seg);

    theta_coord = 2*pi * s ./ total;
    theta_coord = mod(theta_coord - theta_coord(1), 2*pi);

    integrand = f(1:n);
end
