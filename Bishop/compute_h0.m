function h = compute_h0(R, Z, varargin)
%COMPUTE_H0 Compute h0(l) from Bishop Eq. (14).
%
%   h = COMPUTE_H0(R, Z)
%   h = COMPUTE_H0(R, Z, 'Name', value, ...)
%
% Inputs
%   R, Z : coordinates of the reference flux surface, ordered along l
%
% Name-value options
%   'Closed'            : true/false/[], default []
%   'Method'            : 'gradient' or 'central', default 'gradient'
%   'Direction'         : 'as-is', 'clockwise', or 'counterclockwise'
%                         default 'clockwise'
%   'RemoveDuplicateLastPoint' : true/false, default true
%   'X0'                : scalar, optional.
%                         If omitted, use R(1) after preprocessing.
%   'UseIntegralForm'   : true/false, default false
%                         false -> use h0 = R/R0   (recommended)
%                         true  -> use h0 = 1 + (1/R0) * int cos(u) dl
%
% Output fields
%   h.R, h.Z
%   h.l, h.L
%   h.X0
%   h.u
%   h.cosu
%   h.h0
%   h.h0_from_R
%   h.h0_from_integral
%
% Notes
%   In Bishop's notation X is the cylindrical major-radius coordinate.
%   In GEQDSK notation that is usually R.

    p = inputParser;
    addParameter(p, 'Closed', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Method', 'gradient', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Direction', 'clockwise', @(x)ischar(x) || isstring(x));
    addParameter(p, 'RemoveDuplicateLastPoint', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'X0', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'UseIntegralForm', false, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    % geometry and u
    geom = compute_tangent_normal(R, Z, ...
        'Closed', opt.Closed, ...
        'Method', opt.Method, ...
        'NormalDirection', 'bishop', ...
        'Direction', opt.Direction, ...
        'RemoveDuplicateLastPoint', opt.RemoveDuplicateLastPoint);

    R = geom.R(:);
    Z = geom.Z(:);
    l = geom.l(:);

    u = unwrap(atan2(geom.tZ, geom.tR));
    cosu = cos(u);

    if isempty(opt.X0)
        X0 = R(1);
    else
        X0 = opt.X0;
    end

    if abs(X0) < eps
        error('compute_h0:BadX0', 'X0 is too close to zero.');
    end

    % exact discrete version from X(l)/X0
    h0_from_R = R / X0;

    % integral form from Eq. (14)
    h0_from_integral = 1 + cumtrapz(l, cosu) / X0;

    if opt.UseIntegralForm
        h0 = h0_from_integral;
    else
        h0 = h0_from_R;
    end

    h = struct();
    h.R = R;
    h.Z = Z;
    h.l = l;
    h.L = geom.L;
    h.X0 = X0;

    h.u = u;
    h.cosu = cosu;

    h.h0 = h0;
    h.h0_from_R = h0_from_R;
    h.h0_from_integral = h0_from_integral;
end
