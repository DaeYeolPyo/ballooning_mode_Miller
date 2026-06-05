function curv = compute_curvature(R, Z, varargin)
%COMPUTE_CURVATURE Compute signed curvature and curvature radius.
%
%   curv = COMPUTE_CURVATURE(R, Z)
%   curv = COMPUTE_CURVATURE(R, Z, 'Name', value, ...)
%
% Inputs
%   R, Z : vectors of curve coordinates ordered along the curve
%
% Name-value options
%   'Closed'            : true/false/[], default []
%   'Method'            : 'gradient' or 'central', default 'central'
%   'NormalDirection'   : 'left', 'right', 'outward', 'inward', 'bishop'
%                         default 'bishop'
%   'Direction'         : 'as-is', 'clockwise', or 'counterclockwise'
%                         default 'as-is'
%   'RemoveDuplicateLastPoint' : true/false, default true
%
% Output struct fields
%   curv.R, curv.Z
%   curv.l              : cumulative arc length on stored points
%   curv.L              : total closed-curve length including closing segment
%   curv.tR, curv.tZ    : unit tangent
%   curv.nR, curv.nZ    : selected unit normal
%   curv.u              : unwrapped tangent angle
%   curv.dtRdl, dtZdl
%   curv.dudl
%   curv.kappa_signed   : signed curvature, dt/dl = kappa_signed*n
%   curv.kappa          : abs(kappa_signed)
%   curv.Rc_signed      : 1/kappa_signed
%   curv.Rc             : 1/kappa
%
% Bishop convention:
%   With Direction='clockwise' and NormalDirection='bishop',
%   n = (sin u, -cos u) and R(l,rho) = R(l) + rho*sin(u).

    p = inputParser;
    addParameter(p, 'Closed', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Method', 'central', @(x)ischar(x) || isstring(x));
    addParameter(p, 'NormalDirection', 'bishop', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Direction', 'as-is', @(x)ischar(x) || isstring(x));
    addParameter(p, 'RemoveDuplicateLastPoint', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    geom = compute_tangent_normal(R, Z, ...
        'Closed', opt.Closed, ...
        'Method', opt.Method, ...
        'NormalDirection', opt.NormalDirection, ...
        'Direction', opt.Direction, ...
        'RemoveDuplicateLastPoint', opt.RemoveDuplicateLastPoint);

    R = geom.R(:);
    Z = geom.Z(:);
    l = geom.l(:);
    npts = numel(l);

    if geom.closed
        dl = geom.dl(:);
        [dtRdl, dtZdl] = closed_central_derivative(geom.tR, geom.tZ, dl);
        dudl = closed_angle_derivative(geom.u, dl);
    else
        dtRdl = gradient(geom.tR, l);
        dtZdl = gradient(geom.tZ, l);
        dudl = gradient(unwrap(geom.u), l);
    end

    kappa_signed = dtRdl .* geom.nR + dtZdl .* geom.nZ;
    kappa = abs(kappa_signed);

    Rc_signed = inf(npts, 1);
    goodSigned = abs(kappa_signed) > eps & isfinite(kappa_signed);
    Rc_signed(goodSigned) = 1 ./ kappa_signed(goodSigned);

    Rc = inf(npts, 1);
    goodAbs = kappa > eps & isfinite(kappa);
    Rc(goodAbs) = 1 ./ kappa(goodAbs);

    curv = struct();
    curv.R = R;
    curv.Z = Z;
    curv.l = l;
    curv.L = geom.L;
    curv.dl = geom.dl;
    curv.closed = geom.closed;

    curv.tR = geom.tR;
    curv.tZ = geom.tZ;
    curv.nR = geom.nR;
    curv.nZ = geom.nZ;
    curv.u = geom.u;
    curv.dudl = dudl;

    curv.dtRdl = dtRdl;
    curv.dtZdl = dtZdl;

    curv.kappa_signed = kappa_signed;
    curv.kappa = kappa;
    curv.Rc_signed = Rc_signed;
    curv.Rc = Rc;

    curv.direction = geom.direction;
    curv.normalDirection = geom.normalDirection;
end

%======================================================================
function [dfRdl, dfZdl] = closed_central_derivative(fR, fZ, dl)
    fR = fR(:);
    fZ = fZ(:);
    dl = dl(:);

    npts = numel(fR);
    dfRdl = zeros(npts, 1);
    dfZdl = zeros(npts, 1);

    for i = 1:npts
        im = i - 1; if im < 1, im = npts; end
        ip = i + 1; if ip > npts, ip = 1; end

        ds = dl(im) + dl(i);
        dfRdl(i) = (fR(ip) - fR(im)) / ds;
        dfZdl(i) = (fZ(ip) - fZ(im)) / ds;
    end
end

%======================================================================
function dudl = closed_angle_derivative(u, dl)
    u = unwrap(u(:));
    dl = dl(:);

    npts = numel(u);
    dudl = zeros(npts, 1);

    for i = 1:npts
        im = i - 1; if im < 1, im = npts; end
        ip = i + 1; if ip > npts, ip = 1; end

        ui = u;
        if ip == 1
            ui_ip = u(1) + 2*pi*round((u(end) - u(1))/(2*pi));
        else
            ui_ip = ui(ip);
        end

        if im == npts
            ui_im = u(end) - 2*pi*round((u(end) - u(1))/(2*pi));
        else
            ui_im = ui(im);
        end

        ds = dl(im) + dl(i);
        dudl(i) = (ui_ip - ui_im) / ds;
    end
end
