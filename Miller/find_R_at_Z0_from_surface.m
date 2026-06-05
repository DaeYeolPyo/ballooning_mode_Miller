function R0 = find_R_at_Z0_from_surface(surf, varargin)
%FIND_R_AT_Z0_FROM_SURFACE Find R coordinates where extracted flux surface crosses Z=0.
%
%   R0 = find_R_at_Z0_from_surface(surf)
%   R0 = find_R_at_Z0_from_surface(surf, 'Z0', value)
%
% Input
%   surf : struct returned by extract_flux_surface.m
%          must contain surf.R and surf.Z
%
% Options
%   'Z0'        : target Z value, default 0
%   'CloseCurve': true/false, default true
%   'TolZ'      : tolerance for detecting exact Z=Z0 points, default 1e-12
%   'UniqueTol' : tolerance for merging duplicate R crossings, default 1e-8
%   'Sort'      : true/false, default true
%
% Output
%   R0 : R coordinates where the surface crosses Z=Z0.
%
% Notes
%   For a normal closed tokamak flux surface, R0 usually has two values:
%     R0(1) = inboard midplane crossing
%     R0(2) = outboard midplane crossing
%
%   If the surface is shaped strongly, has an X-point-like indentation,
%   or Z0 intersects it more than twice, this function returns all crossings.

    p = inputParser;
    addParameter(p, 'Z0', 0, @(x)isnumeric(x) && isscalar(x));
    addParameter(p, 'CloseCurve', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'TolZ', 1e-12, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'UniqueTol', 1e-8, @(x)isnumeric(x) && isscalar(x) && x >= 0);
    addParameter(p, 'Sort', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    R = surf.R(:);
    Z = surf.Z(:);

    if numel(R) ~= numel(Z)
        error('find_R_at_Z0_fromSurface:SizeMismatch', ...
              'surf.R and surf.Z must have the same number of elements.');
    end

    if numel(R) < 2
        error('find_R_at_Z0_fromSurface:TooFewPoints', ...
              'At least two contour points are required.');
    end

    Zrel = Z - opt.Z0;

    % Close contour if requested and not already closed.
    if opt.CloseCurve
        if hypot(R(end) - R(1), Z(end) - Z(1)) > 1e-14
            R = [R; R(1)];
            Zrel = [Zrel; Zrel(1)];
        end
    end

    Rcross = [];

    for i = 1:numel(R)-1
        R1 = R(i);
        R2 = R(i+1);

        Z1 = Zrel(i);
        Z2 = Zrel(i+1);

        % Case 1: point exactly on Z=Z0
        if abs(Z1) <= opt.TolZ
            Rcross(end+1,1) = R1; %#ok<AGROW>
        end

        % Case 2: segment crosses Z=Z0
        if Z1 * Z2 < 0
            % Linear interpolation along segment
            t = -Z1 / (Z2 - Z1);
            Rc = R1 + t * (R2 - R1);

            Rcross(end+1,1) = Rc; %#ok<AGROW>
        end
    end

    % Check final point if exactly on Z0
    if abs(Zrel(end)) <= opt.TolZ
        Rcross(end+1,1) = R(end); %#ok<AGROW>
    end

    if isempty(Rcross)
        R0 = [];
        return;
    end

    % Remove duplicate crossings caused by exact-zero vertices
    if opt.Sort
        Rcross = sort(Rcross);
    end

    R0 = merge_close_values(Rcross, opt.UniqueTol);
end

% ======================================================================
function xout = merge_close_values(x, tol)
%MERGE_CLOSE_VALUES Merge nearly identical values by averaging clusters.

    x = sort(x(:));

    if isempty(x)
        xout = x;
        return;
    end

    groups = {};
    current = x(1);

    for i = 2:numel(x)
        if abs(x(i) - current(end)) <= tol
            current(end+1,1) = x(i); %#ok<AGROW>
        else
            groups{end+1} = current; %#ok<AGROW>
            current = x(i);
        end
    end

    groups{end+1} = current;

    xout = zeros(numel(groups),1);
    for i = 1:numel(groups)
        xout(i) = mean(groups{i});
    end
end