function geom = compute_tangent_normal(R, Z, varargin)
%COMPUTE_TANGENT_NORMAL Compute unit tangent and normal vectors on a 2D curve.
%
%   geom = COMPUTE_TANGENT_NORMAL(R, Z)
%   geom = COMPUTE_TANGENT_NORMAL(R, Z, 'Name', value, ...)
%
% Inputs
%   R, Z : vectors of curve coordinates ordered along the curve
%
% Name-value options
%   'Closed'            : true/false/[], default []
%                         If [], auto-detect from endpoints.
%   'Method'            : 'gradient' or 'central', default 'gradient'
%                         Method to compute dR/dl and dZ/dl.
%   'NormalDirection'   : 'left', 'right', 'outward', 'inward', 'bishop'
%                         default 'left'
%                         Bishop normal is n = (sin u, -cos u), i.e. the
%                         right normal for t = (cos u, sin u).
%   'Direction'         : 'as-is', 'clockwise', or 'counterclockwise'
%                         default 'as-is'
%   'RemoveDuplicateLastPoint' : true/false, default true
%
% Output struct fields
%   geom.R, geom.Z
%   geom.l              : cumulative arc length
%   geom.L              : total arc length
%   geom.closed
%   geom.tR, geom.tZ    : unit tangent components
%   geom.nR, geom.nZ    : unit normal components
%   geom.u              : tangent angle, atan2(tZ, tR)
%   geom.dRdl, geom.dZdl
%
% Notes
% - 'left' normal means n = (-tZ, tR)
% - 'right' normal means n = ( tZ,-tR)
% - For closed curves, 'outward'/'inward' are inferred from polygon orientation.

    p = inputParser;
    addParameter(p, 'Closed', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'Method', 'gradient', @(x)ischar(x) || isstring(x));
    addParameter(p, 'NormalDirection', 'left', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Direction', 'as-is', @(x)ischar(x) || isstring(x));
    addParameter(p, 'RemoveDuplicateLastPoint', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    %-----------------------------------
    % 1) Arc length preprocessing
    %-----------------------------------
    arc = compute_arclength(R, Z, ...
        'Closed', opt.Closed, ...
        'RemoveDuplicateLastPoint', opt.RemoveDuplicateLastPoint);

    direction = lower(string(opt.Direction));
    if ~any(direction == ["as-is", "clockwise", "counterclockwise"])
        error('compute_tangent_normal:BadDirection', ...
            'Unknown Direction: %s', opt.Direction);
    end

    if arc.closed && direction ~= "as-is"
        Asigned = signed_polygon_area_local(arc.R, arc.Z);
        isClockwise = Asigned < 0;

        if (direction == "clockwise" && ~isClockwise) || ...
           (direction == "counterclockwise" && isClockwise)
            arc = compute_arclength(flipud(arc.R(:)), flipud(arc.Z(:)), ...
                'Closed', true, ...
                'RemoveDuplicateLastPoint', false);
        end
    end

    R = arc.R(:);
    Z = arc.Z(:);
    l = arc.l(:);
    closed = arc.closed;
    npts = arc.npts;

    if npts < 2
        error('compute_tangent_normal:TooFewPoints', ...
            'At least two processed points are required.');
    end

    method = lower(string(opt.Method));

    %-----------------------------------
    % 2) Compute derivatives wrt arc length
    %-----------------------------------
    switch method
        case "gradient"
            if closed
                [dRdl, dZdl] = closed_central_derivative(R, Z, arc.dl);
            else
                dRdl = gradient(R, l);
                dZdl = gradient(Z, l);
            end

        case "central"
            dRdl = zeros(npts,1);
            dZdl = zeros(npts,1);

            if closed
                [dRdl, dZdl] = closed_central_derivative(R, Z, arc.dl);
            else
                % forward/backward at boundaries, central inside
                dRdl(1) = (R(2)-R(1)) / (l(2)-l(1));
                dZdl(1) = (Z(2)-Z(1)) / (l(2)-l(1));

                for i = 2:npts-1
                    ds = l(i+1) - l(i-1);
                    dRdl(i) = (R(i+1)-R(i-1)) / ds;
                    dZdl(i) = (Z(i+1)-Z(i-1)) / ds;
                end

                dRdl(npts) = (R(npts)-R(npts-1)) / (l(npts)-l(npts-1));
                dZdl(npts) = (Z(npts)-Z(npts-1)) / (l(npts)-l(npts-1));
            end

        otherwise
            error('compute_tangent_normal:BadMethod', ...
                'Unknown Method: %s', opt.Method);
    end

    %-----------------------------------
    % 3) Normalize to get unit tangent
    %-----------------------------------
    mag = hypot(dRdl, dZdl);
    bad = mag <= 0 | ~isfinite(mag);

    if any(bad)
        error('compute_tangent_normal:BadDerivative', ...
            'Encountered zero or invalid tangent magnitude.');
    end

    tR = dRdl ./ mag;
    tZ = dZdl ./ mag;

    %-----------------------------------
    % 4) Construct normal
    %-----------------------------------
    ndir = lower(string(opt.NormalDirection));

    switch ndir
        case "left"
            nR = -tZ;
            nZ =  tR;

        case "right"
            nR =  tZ;
            nZ = -tR;

        case "bishop"
            nR =  tZ;
            nZ = -tR;

        case {"outward", "inward"}
            % Start from left normal, then choose sign using polygon orientation
            nR_left = -tZ;
            nZ_left =  tR;

            if ~closed
                warning('compute_tangent_normal:OpenCurveNormal', ...
                    ['NormalDirection="%s" requested for an open curve. ', ...
                     'Using left normal instead.'], ndir);
                nR = nR_left;
                nZ = nZ_left;
            else
                Asigned = signed_polygon_area_local(R, Z);

                % CCW => interior is on left side of tangent
                % so left normal points inward, right normal outward
                if Asigned > 0
                    inwardR  = nR_left;
                    inwardZ  = nZ_left;
                    outwardR = -nR_left;
                    outwardZ = -nZ_left;
                else
                    outwardR = nR_left;
                    outwardZ = nZ_left;
                    inwardR  = -nR_left;
                    inwardZ  = -nZ_left;
                end

                if ndir == "outward"
                    nR = outwardR;
                    nZ = outwardZ;
                else
                    nR = inwardR;
                    nZ = inwardZ;
                end
            end

        otherwise
            error('compute_tangent_normal:BadNormalDirection', ...
                'Unknown NormalDirection: %s', opt.NormalDirection);
    end

    %-----------------------------------
    % 5) Tangent angle
    %-----------------------------------
    u = atan2(tZ, tR);

    %-----------------------------------
    % 6) Output
    %-----------------------------------
    geom = struct();
    geom.R = R;
    geom.Z = Z;
    geom.l = l;
    geom.L = arc.L;
    geom.dl = arc.dl(:);
    geom.closed = closed;

    geom.dRdl = dRdl;
    geom.dZdl = dZdl;

    geom.tR = tR;
    geom.tZ = tZ;

    geom.nR = nR;
    geom.nZ = nZ;

    geom.u = u;
    geom.direction = char(direction);
    geom.normalDirection = char(ndir);
end

%======================================================================
function [dRdl, dZdl] = closed_central_derivative(R, Z, dl)
    R = R(:);
    Z = Z(:);
    dl = dl(:);

    npts = numel(R);
    if numel(dl) ~= npts
        error('compute_tangent_normal:BadClosedDl', ...
            'Closed-curve dl must have the same length as R/Z.');
    end

    dRdl = zeros(npts, 1);
    dZdl = zeros(npts, 1);

    for i = 1:npts
        im = i - 1; if im < 1, im = npts; end
        ip = i + 1; if ip > npts, ip = 1; end

        ds = dl(im) + dl(i);
        dRdl(i) = (R(ip) - R(im)) / ds;
        dZdl(i) = (Z(ip) - Z(im)) / ds;
    end
end

%======================================================================
function A = signed_polygon_area_local(R, Z)
    R = R(:);
    Z = Z(:);

    if R(1) ~= R(end) || Z(1) ~= Z(end)
        R = [R; R(1)];
        Z = [Z; Z(1)];
    end

    A = 0.5 * sum(R(1:end-1).*Z(2:end) - R(2:end).*Z(1:end-1));
end
