function loc = build_local_rho_coordinate(R, Z, rho, varargin)
%BUILD_LOCAL_RHO_COORDINATE Build local (l, rho) coordinates near a flux surface.
%
%   loc = BUILD_LOCAL_RHO_COORDINATE(R, Z, rho)
%   loc = BUILD_LOCAL_RHO_COORDINATE(R, Z, rho, 'Name', value, ...)
%
% Inputs
%   R, Z : reference curve coordinates
%   rho  : scalar or vector of normal offsets [m]
%
% Name-value options
%   'Closed'            : true/false/[], default []
%   'Method'            : 'gradient' or 'central', default 'gradient'
%   'NormalDirection'   : 'left', 'right', 'outward', 'inward', 'bishop'
%                         default 'bishop'
%   'Direction'         : 'as-is', 'clockwise', or 'counterclockwise'
%                         default 'clockwise'
%   'RemoveDuplicateLastPoint' : true/false, default true
%
% Output
%   loc.l               : arc-length coordinate on reference surface
%   loc.rho             : rho values used
%   loc.R0, loc.Z0      : reference curve
%   loc.tR, loc.tZ      : tangent
%   loc.nR, loc.nZ      : normal
%   loc.kappa, loc.Rc   : curvature info
%   loc.R               : size [nl, nrho]
%   loc.Z               : size [nl, nrho]
%
% Convention:
%   R(i,j) = R0(i) + rho(j)*nR(i)
%   Z(i,j) = Z0(i) + rho(j)*nZ(i)
%   With the default Bishop convention, l increases clockwise and
%   n = (sin u, -cos u).

    p = inputParser;
    addParameter(p, 'Closed', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Method', 'gradient', @(x)ischar(x) || isstring(x));
    addParameter(p, 'NormalDirection', 'bishop', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Direction', 'clockwise', @(x)ischar(x) || isstring(x));
    addParameter(p, 'RemoveDuplicateLastPoint', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    rho = rho(:).';  % row vector

    curv = compute_curvature(R, Z, ...
        'Closed', opt.Closed, ...
        'Method', opt.Method, ...
        'NormalDirection', opt.NormalDirection, ...
        'Direction', opt.Direction, ...
        'RemoveDuplicateLastPoint', opt.RemoveDuplicateLastPoint);

    R0 = curv.R(:);
    Z0 = curv.Z(:);
    nR = curv.nR(:);
    nZ = curv.nZ(:);

    nl = numel(R0);
    nrho = numel(rho);

    Rmap = R0 + nR * rho;
    Zmap = Z0 + nZ * rho;

    loc = struct();
    loc.l = curv.l;
    loc.L = curv.L;
    loc.rho = rho;

    loc.R0 = R0;
    loc.Z0 = Z0;

    loc.tR = curv.tR;
    loc.tZ = curv.tZ;
    loc.nR = curv.nR;
    loc.nZ = curv.nZ;

    loc.u = curv.u;
    loc.kappa = curv.kappa;
    loc.kappa_signed = curv.kappa_signed;
    loc.Rc = curv.Rc;
    loc.Rc_signed = curv.Rc_signed;

    loc.R = Rmap;
    loc.Z = Zmap;
    loc.closed = curv.closed;
end
