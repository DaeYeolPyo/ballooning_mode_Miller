function uang = compute_u_angle(R, Z, varargin)
%COMPUTE_U_ANGLE Compute tangent angle u(l) along a flux surface.
%
%   uang = COMPUTE_U_ANGLE(R, Z)
%   uang = COMPUTE_U_ANGLE(R, Z, 'Name', value, ...)
%
% Inputs
%   R, Z : curve coordinates ordered along the flux surface
%
% Name-value options
%   'Closed'            : true/false/[], default []
%   'Method'            : 'gradient' or 'central', default 'gradient'
%   'NormalDirection'   : 'left', 'right', 'outward', 'inward', 'bishop'
%                         default 'bishop'
%   'Direction'         : 'as-is', 'clockwise', or 'counterclockwise'
%                         default 'as-is'
%   'RemoveDuplicateLastPoint' : true/false, default true
%   'Unwrap'            : true/false, default true
%   'ShiftToZero'       : true/false, default false
%                         If true, subtract u(1) so that u(1)=0.
%
% Output fields
%   uang.R, uang.Z
%   uang.l, uang.L
%   uang.tR, uang.tZ
%   uang.nR, uang.nZ
%   uang.u_raw          : atan2(tZ,tR) in [-pi, pi]
%   uang.u              : unwrapped angle (if Unwrap=true)
%   uang.cosu
%   uang.sinu

    p = inputParser;
    addParameter(p, 'Closed', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Method', 'gradient', @(x)ischar(x) || isstring(x));
    addParameter(p, 'NormalDirection', 'bishop', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Direction', 'as-is', @(x)ischar(x) || isstring(x));
    addParameter(p, 'RemoveDuplicateLastPoint', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Unwrap', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'ShiftToZero', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    geom = compute_tangent_normal(R, Z, ...
        'Closed', opt.Closed, ...
        'Method', opt.Method, ...
        'NormalDirection', opt.NormalDirection, ...
        'Direction', opt.Direction, ...
        'RemoveDuplicateLastPoint', opt.RemoveDuplicateLastPoint);

    u_raw = atan2(geom.tZ, geom.tR);

    if opt.Unwrap
        u = unwrap(u_raw);
    else
        u = u_raw;
    end

    if opt.ShiftToZero
        u = u - u(1);
    end

    uang = struct();
    uang.R = geom.R;
    uang.Z = geom.Z;
    uang.l = geom.l;
    uang.L = geom.L;
    uang.closed = geom.closed;
    uang.direction = geom.direction;
    uang.normalDirection = geom.normalDirection;

    uang.tR = geom.tR;
    uang.tZ = geom.tZ;
    uang.nR = geom.nR;
    uang.nZ = geom.nZ;

    uang.u_raw = u_raw;
    uang.u = u;
    uang.cosu = cos(u);
    uang.sinu = sin(u);
end
