function surf = extract_flux_surface(eq, psin_target, varargin)
%EXTRACT_FLUX_SURFACE Extract a single flux surface from GEQDSK equilibrium.
%
%   surf = EXTRACT_FLUX_SURFACE(eq, psin_target)
%   surf = EXTRACT_FLUX_SURFACE(eq, psin_target, 'Name', value, ...)
%
% Inputs
%   eq           : struct from read_geqdsk.m
%   psin_target  : desired normalized poloidal flux in [0,1]
%
% Name-value options
%   'ClosedOnly'     : true/false, default true
%   'MinPoints'      : minimum number of contour points, default 20
%   'SelectMethod'   : 'closest_to_axis' or 'longest', default 'closest_to_axis'
%   'Clockwise'      : true/false/[], default []
%                      If true, reorder points clockwise.
%                      If false, reorder points counterclockwise.
%                      If [], keep as extracted.
%   'UniformTheta'   : true/false, default false
%                      If true, resample selected contour on a uniform theta grid.
%   'NTheta'         : number of theta points, default 256
%   'ThetaCenter'    : 'axis' or 'centroid', default 'axis'
%
% Output struct fields
%   surf.psin
%   surf.psi
%   surf.R
%   surf.Z
%   surf.theta           (when UniformTheta=true)
%   surf.isClosed
%   surf.length
%   surf.area
%   surf.centroidR
%   surf.centroidZ
%   surf.index
%   surf.allContours
%
% Notes
% - Uses contourc on psi(R,Z).
% - Assumes normalized flux is defined as
%       psin = (psi - simag)/(sibry - simag)
%   consistent with standard GEQDSK convention.

    p = inputParser;
    addParameter(p, 'ClosedOnly', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'MinPoints', 20, @(x)isnumeric(x) && isscalar(x) && x >= 2);
    addParameter(p, 'SelectMethod', 'closest_to_axis', @(x)ischar(x) || isstring(x));
    addParameter(p, 'Clockwise', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'UniformTheta', false, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'NTheta', 256, @(x)isnumeric(x) && isscalar(x) && x >= 8);
    addParameter(p, 'ThetaCenter', 'axis', @(x)ischar(x) || isstring(x));
    parse(p, varargin{:});
    opt = p.Results;

    if psin_target < 0 || psin_target > 1
        warning('extract_flux_surface:psinOutOfRange', ...
            'psin_target = %.6g is outside [0,1]. Continuing anyway.', psin_target);
    end

    %-----------------------------------
    % 1) Build psi target from psin
    %-----------------------------------
    dpsi = eq.sibry - eq.simag;
    psi_target = eq.simag + psin_target * dpsi;

    %-----------------------------------
    % 2) Extract all contours at this level
    %-----------------------------------
    C = contourc(eq.rgrid, eq.zgrid, eq.psirz.', [psi_target psi_target]);

    contours = parse_contourc_output(C);

    if isempty(contours)
        error('extract_flux_surface:NoContourFound', ...
            'No contour found for psin_target = %.6f', psin_target);
    end

    %-----------------------------------
    % 3) Compute metrics for each contour
    %-----------------------------------
    keep = false(size(contours));
    for i = 1:numel(contours)
        R = contours(i).R(:);
        Z = contours(i).Z(:);

        npts = numel(R);
        isClosed = is_closed_curve(R, Z);

        len = polyline_length(R, Z);
        [areaVal, cR, cZ] = polygon_metrics(R, Z);

        contours(i).npts = npts;
        contours(i).isClosed = isClosed;
        contours(i).length = len;
        contours(i).area = areaVal;
        contours(i).centroidR = cR;
        contours(i).centroidZ = cZ;

        ok = npts >= opt.MinPoints;
        if opt.ClosedOnly
            ok = ok && isClosed;
        end
        keep(i) = ok;
    end

    contours = contours(keep);

    if isempty(contours)
        error('extract_flux_surface:NoValidContour', ...
            'Contours were found, but none satisfied the selection criteria.');
    end

    %-----------------------------------
    % 4) Select one contour
    %-----------------------------------
    method = lower(string(opt.SelectMethod));

    switch method
        case "closest_to_axis"
            dist2 = arrayfun(@(s) (s.centroidR - eq.rmaxis)^2 + (s.centroidZ - eq.zmaxis)^2, contours);
            [~, idx] = min(dist2);

        case "longest"
            lens = arrayfun(@(s) s.length, contours);
            [~, idx] = max(lens);

        otherwise
            error('extract_flux_surface:BadSelectMethod', ...
                'Unknown SelectMethod: %s', opt.SelectMethod);
    end

    chosen = contours(idx);

    %-----------------------------------
    % 5) Optionally enforce orientation
    %-----------------------------------
    if ~isempty(opt.Clockwise)
        signedA = signed_polygon_area(chosen.R, chosen.Z);
        isCW = signedA < 0;

        if opt.Clockwise && ~isCW
            chosen.R = flipud(chosen.R(:));
            chosen.Z = flipud(chosen.Z(:));
        elseif ~opt.Clockwise && isCW
            chosen.R = flipud(chosen.R(:));
            chosen.Z = flipud(chosen.Z(:));
        end
    end

    %-----------------------------------
    % 6) Optional: resample to uniform theta grid
    %-----------------------------------
    theta_uniform = [];

    if opt.UniformTheta
        centerMode = lower(string(opt.ThetaCenter));
        switch centerMode
            case "axis"
                Rc = eq.rmaxis;
                Zc = eq.zmaxis;
            case "centroid"
                Rc = chosen.centroidR;
                Zc = chosen.centroidZ;
            otherwise
                error('extract_flux_surface:BadThetaCenter', ...
                    'Unknown ThetaCenter: %s', opt.ThetaCenter);
        end

        [Rr, Zr, theta_uniform] = resample_contour_uniform_theta( ...
            chosen.R(:), chosen.Z(:), Rc, Zc, opt.NTheta);

        chosen.R = Rr;
        chosen.Z = Zr;

        if ~isempty(opt.Clockwise)
            signedA = signed_polygon_area(chosen.R, chosen.Z);
            isCW = signedA < 0;

            if opt.Clockwise && ~isCW
                chosen.R = flipud(chosen.R(:));
                chosen.Z = flipud(chosen.Z(:));
                theta_uniform = flipud(theta_uniform(:));
            elseif ~opt.Clockwise && isCW
                chosen.R = flipud(chosen.R(:));
                chosen.Z = flipud(chosen.Z(:));
                theta_uniform = flipud(theta_uniform(:));
            end
        end

        % resampled contour metric refresh
        chosen.isClosed = is_closed_curve(chosen.R, chosen.Z);
        chosen.length   = polyline_length([chosen.R; chosen.R(1)], [chosen.Z; chosen.Z(1)]);
        [chosen.area, chosen.centroidR, chosen.centroidZ] = polygon_metrics(chosen.R, chosen.Z);
    end

    %-----------------------------------
    % 7) Output
    %-----------------------------------
    surf = struct();
    surf.psin = psin_target;
    surf.psi  = psi_target;
    surf.R    = chosen.R(:);
    surf.Z    = chosen.Z(:);
    surf.theta = theta_uniform(:);
    surf.isClosed = chosen.isClosed;
    surf.length   = chosen.length;
    surf.area     = chosen.area;
    surf.centroidR = chosen.centroidR;
    surf.centroidZ = chosen.centroidZ;
    surf.index = idx;
    surf.allContours = contours;
end

%======================================================================
function contours = parse_contourc_output(C)
    contours = struct('level', {}, 'R', {}, 'Z', {});
    k = 1;
    ncol = size(C, 2);

    while k < ncol
        level = C(1, k);
        npts  = C(2, k);

        cols = (k+1):(k+npts);
        if cols(end) > ncol
            break;
        end

        R = C(1, cols).';
        Z = C(2, cols).';

        s = struct();
        s.level = level;
        s.R = R;
        s.Z = Z;

        contours(end+1) = s; %#ok<AGROW>
        k = k + npts + 1;
    end
end

%======================================================================
function tf = is_closed_curve(R, Z)
    if numel(R) < 3
        tf = false;
        return;
    end

    scale = max([max(R)-min(R), max(Z)-min(Z), 1]);
    tol = 1e-6 * scale;

    tf = hypot(R(end)-R(1), Z(end)-Z(1)) < tol;
end

%======================================================================
function L = polyline_length(R, Z)
    dR = diff(R);
    dZ = diff(Z);
    L = sum(hypot(dR, dZ));
end

%======================================================================
function [A, Cx, Cy] = polygon_metrics(R, Z)
    R = R(:);
    Z = Z(:);

    if numel(R) < 3
        A = 0;
        Cx = mean(R, 'omitnan');
        Cy = mean(Z, 'omitnan');
        return;
    end

    if R(1) ~= R(end) || Z(1) ~= Z(end)
        R = [R; R(1)];
        Z = [Z; Z(1)];
    end

    cross = R(1:end-1).*Z(2:end) - R(2:end).*Z(1:end-1);
    Asigned = 0.5 * sum(cross);
    A = abs(Asigned);

    if abs(Asigned) < eps
        Cx = mean(R(1:end-1), 'omitnan');
        Cy = mean(Z(1:end-1), 'omitnan');
    else
        Cx = sum((R(1:end-1)+R(2:end)).*cross) / (6*Asigned);
        Cy = sum((Z(1:end-1)+Z(2:end)).*cross) / (6*Asigned);
    end
end

%======================================================================
function A = signed_polygon_area(R, Z)
    R = R(:);
    Z = Z(:);

    if R(1) ~= R(end) || Z(1) ~= Z(end)
        R = [R; R(1)];
        Z = [Z; Z(1)];
    end

    A = 0.5 * sum(R(1:end-1).*Z(2:end) - R(2:end).*Z(1:end-1));
end

%======================================================================
function [Rout, Zout, theta_uniform] = resample_contour_uniform_theta(R, Z, Rc, Zc, ntheta)
% Resample a closed contour onto a uniform theta grid around (Rc, Zc).

    R = R(:);
    Z = Z(:);

    % remove duplicate endpoint if present
    if numel(R) >= 2 && hypot(R(end)-R(1), Z(end)-Z(1)) < 1e-12
        R = R(1:end-1);
        Z = Z(1:end-1);
    end

    theta = atan2(Z - Zc, R - Rc);
    theta = unwrap(theta);

    % enforce monotonic theta by sorting
    [theta_sorted, idx] = sort(theta);
    R_sorted = R(idx);
    Z_sorted = Z(idx);

    % periodic extension
    theta_ext = [theta_sorted - 2*pi; theta_sorted; theta_sorted + 2*pi];
    R_ext = [R_sorted; R_sorted; R_sorted];
    Z_ext = [Z_sorted; Z_sorted; Z_sorted];

    theta_min = theta_sorted(1);
    theta_uniform = linspace(theta_min, theta_min + 2*pi, ntheta + 1).';
    theta_uniform(end) = [];  % drop duplicate 2pi endpoint

    Rout = interp1(theta_ext, R_ext, theta_uniform, 'spline');
    Zout = interp1(theta_ext, Z_ext, theta_uniform, 'spline');
end
