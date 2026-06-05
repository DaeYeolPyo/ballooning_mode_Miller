function out = resample_closed_curve_arclength(R, Z, npts, varargin)
%RESAMPLE_CLOSED_CURVE_ARCLENGTH Resample a closed curve on uniform l.
%
%   out = RESAMPLE_CLOSED_CURVE_ARCLENGTH(R, Z, npts)
%   out = RESAMPLE_CLOSED_CURVE_ARCLENGTH(..., 'Name', value, ...)
%
% Name-value options
%   'Direction'    : 'as-is', 'clockwise', or 'counterclockwise'
%                    default 'clockwise'
%   'InterpMethod' : 'pchip', 'spline', or 'linear', default 'pchip'

    p = inputParser;
    addParameter(p, 'Direction', 'clockwise', @(x)ischar(x) || isstring(x));
    addParameter(p, 'InterpMethod', 'pchip', @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    opt = p.Results;

    if npts < 8
        error('resample_closed_curve_arclength:TooFewPoints', ...
            'npts must be at least 8.');
    end

    geom = compute_tangent_normal(R, Z, ...
        'Closed', true, ...
        'Direction', opt.Direction, ...
        'NormalDirection', 'bishop', ...
        'RemoveDuplicateLastPoint', true);

    lExt = [geom.l(:); geom.L];
    RExt = [geom.R(:); geom.R(1)];
    ZExt = [geom.Z(:); geom.Z(1)];

    lUniform = linspace(0, geom.L, round(npts) + 1).';
    lUniform(end) = [];

    method = char(lower(string(opt.InterpMethod)));
    Rout = interp1(lExt, RExt, lUniform, method);
    Zout = interp1(lExt, ZExt, lUniform, method);

    out = struct();
    out.R = Rout(:);
    out.Z = Zout(:);
    out.l = lUniform(:);
    out.L = geom.L;
    out.theta_ballooning = 2*pi * out.l / out.L;
    out.direction = char(lower(string(opt.Direction)));
end
