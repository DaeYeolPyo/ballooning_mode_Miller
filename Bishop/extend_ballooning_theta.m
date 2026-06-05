function ext = extend_ballooning_theta(base, varargin)
%EXTEND_BALLOONING_THETA Periodically extend local geometry in ballooning theta.
%
%   ext = EXTEND_BALLOONING_THETA(base)
%   ext = EXTEND_BALLOONING_THETA(base, 'Name', value, ...)
%
% The ballooning representation treats the poloidal coordinate as an
% extended coordinate theta in (-inf, inf). Numerically this is represented
% by repeating one closed surface over a finite window.
%
% Required base fields
%   base.l : arclength coordinate on one period
%   base.L : total period length
%
% Name-value options
%   'NPeriodsEachSide' : integer, default 3. The output covers
%                        approximately [-2*pi*N, 2*pi*N).
%   'Direction'        : 'auto', 'clockwise', or 'counterclockwise'.
%                        default 'auto'
%   'ThetaField'       : field to use if present, default 'theta_ballooning'
%
% Output
%   ext.theta          : extended theta coordinate
%   ext.l_extended     : extended arclength coordinate
%   ext.turn           : integer period index for each point
%
% Numeric vector fields with the same length as base.l are repeated. The
% unwrapped tangent angle base.u is shifted by +/-2*pi per period so that it
% stays continuous with the chosen l direction.

    p = inputParser;
    addParameter(p, 'NPeriodsEachSide', 3, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'Direction', 'auto', @(x)ischar(x) || isstring(x));
    addParameter(p, 'ThetaField', 'theta_ballooning', @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    opt = p.Results;

    lBase = base.l(:);
    n = numel(lBase);

    if ~isfield(base, 'L') || ~isscalar(base.L) || base.L <= 0
        error('extend_ballooning_theta:BadPeriod', ...
            'base must contain a positive scalar L.');
    end

    thetaField = char(opt.ThetaField);
    if isfield(base, thetaField) && numel(base.(thetaField)) == n
        thetaValue = base.(thetaField);
        thetaBase = thetaValue(:);
    else
        thetaBase = 2*pi * (lBase - lBase(1)) ./ base.L;
    end

    direction = lower(string(opt.Direction));
    if direction == "auto"
        if isfield(base, 'direction')
            direction = lower(string(base.direction));
        else
            direction = infer_direction_from_u(base);
        end
    end

    if ~any(direction == ["clockwise", "counterclockwise", "as-is"])
        error('extend_ballooning_theta:BadDirection', ...
            'Unknown Direction: %s', opt.Direction);
    end

    if direction == "counterclockwise"
        uShiftSign = 1;
    else
        % Bishop report convention: l increases clockwise, so unwrapped u
        % decreases by 2*pi per poloidal turn.
        uShiftSign = -1;
    end

    nside = round(opt.NPeriodsEachSide);
    turns = (-nside):(nside-1);

    theta = zeros(n * numel(turns), 1);
    lExtended = zeros(size(theta));
    turnIndex = zeros(size(theta));

    cursor = 1;
    for k = turns
        idx = cursor:(cursor+n-1);
        theta(idx) = thetaBase + 2*pi*k;
        lExtended(idx) = lBase + base.L*k;
        turnIndex(idx) = k;
        cursor = cursor + n;
    end

    [theta, ord] = sort(theta);
    lExtended = lExtended(ord);
    turnIndex = turnIndex(ord);

    ext = struct();
    ext.theta = theta;
    ext.l_extended = lExtended;
    ext.turn = turnIndex;
    ext.L = base.L;
    ext.direction = char(direction);
    ext.NPeriodsEachSide = nside;

    names = fieldnames(base);
    for i = 1:numel(names)
        name = names{i};
        value = base.(name);

        if strcmp(name, 'l')
            ext.l = lExtended;
            continue;
        end

        if strcmp(name, 'L')
            continue;
        end

        if isnumeric(value) && isvector(value) && numel(value) == n
            repeated = repeat_vector(value(:), turns, name, uShiftSign);
            ext.(name) = repeated(ord);
        elseif isnumeric(value) && size(value, 1) == n
            repeated = repeat_rows(value, turns, name, uShiftSign);
            ext.(name) = repeated(ord, :);
        else
            ext.(name) = value;
        end
    end

    if isfield(ext, 'u')
        ext.cosu = cos(ext.u);
        ext.sinu = sin(ext.u);
    end
end

%======================================================================
function repeated = repeat_vector(value, turns, name, uShiftSign)
    n = numel(value);
    repeated = zeros(n * numel(turns), 1);

    cursor = 1;
    for k = turns
        idx = cursor:(cursor+n-1);
        if strcmp(name, 'u')
            repeated(idx) = value + uShiftSign * 2*pi*k;
        else
            repeated(idx) = value;
        end
        cursor = cursor + n;
    end
end

%======================================================================
function repeated = repeat_rows(value, turns, name, uShiftSign)
    [n, m] = size(value);
    repeated = zeros(n * numel(turns), m);

    cursor = 1;
    for k = turns
        idx = cursor:(cursor+n-1);
        if strcmp(name, 'u')
            repeated(idx, :) = value + uShiftSign * 2*pi*k;
        else
            repeated(idx, :) = value;
        end
        cursor = cursor + n;
    end
end

%======================================================================
function direction = infer_direction_from_u(base)
    direction = "clockwise";

    if ~isfield(base, 'u')
        return;
    end

    u = unwrap(base.u(:));
    if numel(u) < 2
        return;
    end

    if u(end) > u(1)
        direction = "counterclockwise";
    end
end
