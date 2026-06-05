function arc = compute_arclength(R, Z, varargin)
%COMPUTE_ARCLENGTH Compute cumulative arc length along a 2D curve.
%
%   arc = COMPUTE_ARCLENGTH(R, Z)
%   arc = COMPUTE_ARCLENGTH(R, Z, 'Name', value, ...)
%
% Inputs
%   R, Z : vectors of curve coordinates
%
% Name-value options
%   'Closed'        : true/false/[], default []
%                     If [], auto-detect from endpoints.
%                     If true, enforce closed-curve handling.
%                     If false, treat as open curve.
%
%   'CloseTolerance': endpoint tolerance for closure detection,
%                     default = 1e-6 * characteristic scale
%
%   'RemoveDuplicateLastPoint' : true/false, default true
%                     If curve is closed and last point duplicates first,
%                     remove the last point before computing arc length.
%
% Outputs (struct)
%   arc.R           : processed R points
%   arc.Z           : processed Z points
%   arc.npts        : number of processed points
%   arc.closed      : whether treated as closed
%   arc.dl          : segment lengths. For a closed curve this includes
%                     the final segment from the last point back to the
%                     first point.
%   arc.l           : cumulative arc length, same length as R/Z
%   arc.L           : total length
%
% Example
%   surf = extract_flux_surface(eq, 0.5);
%   arc  = compute_arclength(surf.R, surf.Z);
%   plot(arc.l, arc.R)

    p = inputParser;
    addParameter(p, 'Closed', [], @(x)islogical(x) && isscalar(x) || isempty(x));
    addParameter(p, 'CloseTolerance', [], @(x)isnumeric(x) && isscalar(x) || isempty(x));
    addParameter(p, 'RemoveDuplicateLastPoint', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    R = R(:);
    Z = Z(:);

    if numel(R) ~= numel(Z)
        error('compute_arclength:SizeMismatch', ...
            'R and Z must have the same number of elements.');
    end

    if numel(R) < 2
        error('compute_arclength:TooFewPoints', ...
            'At least two points are required.');
    end

    %-----------------------------------
    % 1) Determine closure
    %-----------------------------------
    scale = max([max(R)-min(R), max(Z)-min(Z), 1]);
    if isempty(opt.CloseTolerance)
        tol = 1e-6 * scale;
    else
        tol = opt.CloseTolerance;
    end

    autoClosed = hypot(R(end) - R(1), Z(end) - Z(1)) < tol;

    if isempty(opt.Closed)
        isClosed = autoClosed;
    else
        isClosed = opt.Closed;
    end

    %-----------------------------------
    % 2) If closed and duplicated endpoint exists, remove it
    %-----------------------------------
    if isClosed && opt.RemoveDuplicateLastPoint
        if hypot(R(end) - R(1), Z(end) - Z(1)) < tol
            R = R(1:end-1);
            Z = Z(1:end-1);
        end
    end

    npts = numel(R);

    if npts < 2
        error('compute_arclength:TooFewProcessedPoints', ...
            'Too few points remain after processing.');
    end

    %-----------------------------------
    % 3) Segment lengths and cumulative arc length
    %-----------------------------------
    if isClosed
        dR = [diff(R); R(1) - R(end)];
        dZ = [diff(Z); Z(1) - Z(end)];
        dl = hypot(dR, dZ);

        % l is attached to the stored unique points only. The closing
        % segment contributes to L but does not create a duplicate endpoint.
        l = [0; cumsum(dl(1:end-1))];
        L = sum(dl);
    else
        dR = diff(R);
        dZ = diff(Z);
        dl = hypot(dR, dZ);

        l = [0; cumsum(dl)];
        L = l(end);
    end

    %-----------------------------------
    % 5) Output
    %-----------------------------------
    arc = struct();
    arc.R = R;
    arc.Z = Z;
    arc.npts = npts;
    arc.closed = isClosed;
    arc.dl = dl;
    arc.l = l;
    arc.L = L;
end
