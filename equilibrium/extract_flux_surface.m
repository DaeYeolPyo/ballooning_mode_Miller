function surf = extract_flux_surface(eq, eqfunc, target_psiN, opts)
    arguments
        eq struct
        eqfunc struct
        target_psiN (1, 1) double
        opts.Npoints (1, 1) int32 = 1024
        opts.Start string = "outboard"
        opts.Orientation string = "outboard_up"
    end

    Rgrid = eq.rgrid(:);
    Zgrid = eq.zgrid(:);

    target_psi = eq.simag + target_psiN*(eq.sibry - eq.simag);
    C = contourc(Rgrid, Zgrid, eq.psirz.', [target_psi target_psi]);

    segments = parse_contourc_segments(C);

    if isempty(segments)
        error('No contour found for target_psiN = %.16g.', target_psiN);
    end

    idx = choose_axis_enclosing_segment(segments, eq.rmaxis, eq.zmaxis);

    if isempty(idx)
        error('No closed contour enclosing magnetic axis was found.');
    end

    R = segments(idx).R(:);
    Z = segments(idx).Z(:);

    [R, Z] = remove_duplicate_endpoint(R, Z);
    [R, Z] = orient_surface(R, Z, opts.Orientation);
    [R, Z] = start_surface_at_outboard(R, Z, eq.rmaxis, eq.zmaxis);

    [R, Z, s, dl, L] = resample_closed_curve(R, Z, opts.Npoints);

    surf = struct();
    surf.psiN = target_psiN;
    surf.psi = target_psiN*(eq.sibry - eq.simag) + eq.simag;

    surf.R = R;
    surf.Z = Z;
    surf.area = s;
    surf.dl = dl;
    surf.L = L;

    surf.F = eqfunc.F(target_psiN);
    surf.q = eqfunc.q(target_psiN);
    surf.pprime = eqfunc.pprime(target_psiN);
    surf.qprime = eqfunc.qprime(target_psiN);

    surf.Raxis = eq.rmaxis;
    surf.Zaxis = eq.zmaxis;
end

function segments = parse_contourc_segments(C)
    segments = struct('level', {}, 'R', {}, 'Z', {}, 'area', {}, 'closed', {});

    k = 1;
    nseg = 0;

    while k < size(C, 2)
        level = C(1, k);
        n = C(2, k);

        cols = k + (1:n);
        R = C(1, cols).';
        Z = C(2, cols).';

        [R, Z] = remove_duplicate_endpoint(R, Z);

        if numel(R) >= 4
            closed = hypot(R(1) - R(end), Z(1) - Z(end)) < ...
                     1e-8 * max(1, max(hypot(R, Z)));

            % contourc usually gives closed contours with duplicated endpoint,
            % which was removed above. For selection, polygon closure is implied.
            if ~closed
                closed = true;
            end

            nseg = nseg + 1;
            segments(nseg).level = level;
            segments(nseg).R = R;
            segments(nseg).Z = Z;
            segments(nseg).area = polyarea(R, Z);
            segments(nseg).closed = closed;
        end

        k = k + n + 1;
    end
end

function idx = choose_axis_enclosing_segment(segments, Raxis, Zaxis)
    idx = [];

    bestArea = -inf;

    for k = 1:numel(segments)
        R = segments(k).R;
        Z = segments(k).Z;

        if ~segments(k).closed
            continue
        end

        inside = inpolygon(Raxis, Zaxis, R, Z);

        if inside && segments(k).area > bestArea
            bestArea = segments(k).area;
            idx = k;
        end
    end
end

function [R, Z] = remove_duplicate_endpoint(R, Z)
    if numel(R) > 1 && hypot(R(1) - R(end), Z(1) - Z(end)) < ...
            1e-10 * max(1, max(hypot(R, Z)))
        R(end) = [];
        Z(end) = [];
    end
end

function [R, Z] = orient_surface(R, Z, orientation)
    signedArea = signed_polygon_area(R, Z);

    % Positive signed area means counter-clockwise in the R-Z plane.
    isCCW = signedArea > 0;

    switch orientation
        case "ccw"
            if ~isCCW
                R = flipud(R);
                Z = flipud(Z);
            end

        case "cw"
            if isCCW
                R = flipud(R);
                Z = flipud(Z);
            end

        case "outboard_up"
            % First force CCW. After starting at outboard midplane,
            % this usually moves upward first for standard R-horizontal, Z-vertical plots.
            if ~isCCW
                R = flipud(R);
                Z = flipud(Z);
            end

        otherwise
            error('extract_flux_surface:BadOrientation', ...
                  'Unknown orientation option: %s', orientation);
    end
end

function A = signed_polygon_area(R, Z)
    Rp = [R(2:end); R(1)];
    Zp = [Z(2:end); Z(1)];
    A = 0.5 * sum(R .* Zp - Rp .* Z);
end

function [R, Z] = start_surface_at_outboard(R, Z, Raxis, Zaxis)
    % Prefer the point closest to the outboard midplane:
    % large R and Z close to magnetic-axis height.
    Rspan = max(R) - min(R);
    Zspan = max(Z) - min(Z);

    score = ((R - max(R)) / max(Rspan, eps)).^2 + ...
            ((Z - Zaxis) / max(Zspan, eps)).^2;

    [~, i0] = min(score);

    R = circshift(R, 1 - i0);
    Z = circshift(Z, 1 - i0);

    % Make the first step go upward from the outboard side.
    if numel(Z) >= 2 && Z(2) < Z(1)
        R = [R(1); flipud(R(2:end))];
        Z = [Z(1); flipud(Z(2:end))];
    end
end

function [Rq, Zq, sq, dlq, L] = resample_closed_curve(R, Z, N)
    Rclosed = [R; R(1)];
    Zclosed = [Z; Z(1)];

    ds = hypot(diff(Rclosed), diff(Zclosed));
    s = [0; cumsum(ds)];
    L = s(end);

    keep = [true; diff(s) > 0];
    s = s(keep);
    Rclosed = Rclosed(keep);
    Zclosed = Zclosed(keep);

    sq = linspace(0, L, N + 1).';
    sq(end) = [];

    Rq = interp1(s, Rclosed, sq, 'pchip');
    Zq = interp1(s, Zclosed, sq, 'pchip');

    Rq_next = circshift(Rq, -1);
    Zq_next = circshift(Zq, -1);

    dlq = hypot(Rq_next - Rq, Zq_next - Zq);
end