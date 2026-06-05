function surf = build_local_flux_surface_from_geqdsk(filename, psin_target, opts)
%BUILD_LOCAL_FLUX_SURFACE_FROM_GEQDSK
%
% Construct Bishop-style local flux-surface coordinates (l,rho)
% from a GEQDSK equilibrium.
%
% Required:
%   eq = read_geqdsk(filename)
%
% Usage:
%   opts = struct();
%   opts.N = 1024;
%   opts.smoothCutoff = 50;
%   opts.rho = linspace(-0.03, 0.03, 41);
%   opts.plot = true;
%
%   surf = build_local_flux_surface_from_geqdsk('gfile', 0.8, opts);
%
% Outputs:
%   surf.l          : arclength coordinate, 0 <= l < L
%   surf.X          : major-radius coordinate X(l) = R(l) in GEQDSK sense
%   surf.Z          : vertical coordinate Z(l)
%   surf.u          : tangent angle, dX/dl = cos u, dZ/dl = sin u
%   surf.tX, surf.tZ: unit tangent components
%   surf.nX, surf.nZ: Bishop normal = (sin u, -cos u)
%   surf.kappa      : signed curvature
%   surf.Rc_signed  : signed curvature radius in Bishop convention
%   surf.Rc_abs     : absolute curvature radius
%   surf.Xrho,Zrho  : local coordinate mesh X(l,rho), Z(l,rho)
%
% Notes:
%   In Bishop's notation the cylindrical major radius is often written x.
%   Here I use X for that coordinate to avoid confusing it with curvature
%   radius Rc.

    if nargin < 3
        opts = struct();
    end

    opts = set_default_opts(opts);

    % ------------------------------------------------------------
    % 1. Read GEQDSK
    % ------------------------------------------------------------
    eq = read_geqdsk(filename);

    Rgrid = eq.rgrid(:).';
    Zgrid = eq.zgrid(:);

    % GEQDSK reader stores psirz as [nw, nh].
    % contourc expects matrix size [length(Zgrid), length(Rgrid)].
    psiRZ = eq.psirz.';

    psi_target = eq.simag + psin_target * (eq.sibry - eq.simag);

    % ------------------------------------------------------------
    % 2. Extract the selected flux contour
    % ------------------------------------------------------------
    contours = extract_contours(Rgrid, Zgrid, psiRZ, psi_target);

    if isempty(contours)
        error('No contour found for psin_target = %.6g', psin_target);
    end

    [Xraw, Zraw] = choose_flux_surface(contours, eq);

    % Remove duplicate/near-duplicate points.
    [Xraw, Zraw] = clean_closed_curve(Xraw, Zraw);

    % Enforce Bishop's clockwise l direction.
    % With this direction the Bishop normal n = (sin u, -cos u) is the
    % right normal of the curve.
    if signed_polygon_area(Xraw, Zraw) > 0
        Xraw = flipud(Xraw(:));
        Zraw = flipud(Zraw(:));
    else
        Xraw = Xraw(:);
        Zraw = Zraw(:);
    end

    % ------------------------------------------------------------
    % 3. Uniform arclength resampling
    % ------------------------------------------------------------
    [l0, Xu, Zu, Ltot] = resample_closed_curve_by_arclength( ...
        Xraw, Zraw, opts.N);

    % ------------------------------------------------------------
    % 4. Periodic Fourier smoothing and derivatives
    % ------------------------------------------------------------
    [X, dXdl, d2Xdl2] = periodic_fourier_smooth_derivatives( ...
        Xu, Ltot, opts.smoothCutoff, opts.filterPower);

    [Z, dZdl, d2Zdl2] = periodic_fourier_smooth_derivatives( ...
        Zu, Ltot, opts.smoothCutoff, opts.filterPower);

    % Re-normalize tangent because smoothing may slightly perturb arclength.
    speed = sqrt(dXdl.^2 + dZdl.^2);
    tX = dXdl ./ speed;
    tZ = dZdl ./ speed;

    % Tangent angle u:
    %   dX/dl = cos u
    %   dZ/dl = sin u
    u = unwrap(atan2(tZ, tX));

    % Curvature:
    %   kappa = d u / d l
    % equivalently:
    %   kappa = (X' Z'' - Z' X'') / (X'^2 + Z'^2)^(3/2)
    kappa = (dXdl .* d2Zdl2 - dZdl .* d2Xdl2) ./ ...
            (dXdl.^2 + dZdl.^2).^(3/2);

    % Bishop's coordinate uses normal n = (sin u, -cos u).
    % For clockwise l, the geometric signed curvature du/dl is negative
    % for a circular surface, and h_l = 1 - rho/Rc corresponds to
    % Rc_signed = -1/kappa.
    Rc_signed = -1 ./ kappa;
    Rc_abs    = abs(1 ./ kappa);

    % Bishop normal direction
    nX = sin(u);
    nZ = -cos(u);

    % Optional rho mesh
    rho = opts.rho(:).';
    Xrho = X(:) + nX(:) * rho;
    Zrho = Z(:) + nZ(:) * rho;

    % ------------------------------------------------------------
    % 5. Pack output
    % ------------------------------------------------------------
    surf = struct();
    surf.eq = eq;
    surf.psin_target = psin_target;
    surf.psi_target = psi_target;

    surf.l = l0(:);
    surf.L = Ltot;

    surf.Xraw = Xraw(:);
    surf.Zraw = Zraw(:);

    surf.X = X(:);
    surf.Z = Z(:);

    surf.u = u(:);
    surf.tX = tX(:);
    surf.tZ = tZ(:);
    surf.nX = nX(:);
    surf.nZ = nZ(:);

    surf.kappa = kappa(:);
    surf.Rc_signed = Rc_signed(:);
    surf.Rc_abs = Rc_abs(:);

    surf.rho = rho;
    surf.Xrho = Xrho;
    surf.Zrho = Zrho;

    surf.opts = opts;

    if opts.plot
        plot_local_surface_result(surf);
    end
end

% ======================================================================
function opts = set_default_opts(opts)
    if ~isfield(opts, 'N')
        opts.N = 1024;
    end

    if ~isfield(opts, 'smoothCutoff')
        % Number of Fourier modes effectively retained.
        % Increase this if the flux surface is strongly shaped.
        % Decrease this if curvature Rc(l) is still noisy.
        opts.smoothCutoff = 50;
    end

    if ~isfield(opts, 'filterPower')
        opts.filterPower = 8;
    end

    if ~isfield(opts, 'rho')
        opts.rho = linspace(-0.02, 0.02, 31);
    end

    if ~isfield(opts, 'plot')
        opts.plot = true;
    end
end

% ======================================================================
function contours = extract_contours(Rgrid, Zgrid, psiRZ, psi_target)
    C = contourc(Rgrid, Zgrid, psiRZ, [psi_target psi_target]);

    contours = {};
    idx = 1;

    while idx < size(C, 2)
        level = C(1, idx); %#ok<NASGU>
        npts  = C(2, idx);

        pts = C(:, idx+1:idx+npts);
        R = pts(1, :).';
        Z = pts(2, :).';

        contours{end+1} = struct('R', R, 'Z', Z); %#ok<AGROW>

        idx = idx + npts + 1;
    end
end

% ======================================================================
function [X, Z] = choose_flux_surface(contours, eq)
    % Prefer a closed contour containing the magnetic axis.
    bestIdx = [];
    bestArea = -inf;

    for i = 1:numel(contours)
        R = contours{i}.R(:);
        Z = contours{i}.Z(:);

        if numel(R) < 20
            continue;
        end

        % Close if necessary for inpolygon/area checks
        if hypot(R(1)-R(end), Z(1)-Z(end)) > 1e-10
            Rtest = [R; R(1)];
            Ztest = [Z; Z(1)];
        else
            Rtest = R;
            Ztest = Z;
        end

        containsAxis = inpolygon(eq.rmaxis, eq.zmaxis, Rtest, Ztest);
        areaAbs = abs(signed_polygon_area(Rtest, Ztest));

        if containsAxis && areaAbs > bestArea
            bestArea = areaAbs;
            bestIdx = i;
        end
    end

    % Fallback: largest closed-ish contour
    if isempty(bestIdx)
        for i = 1:numel(contours)
            R = contours{i}.R(:);
            Z = contours{i}.Z(:);
            areaAbs = abs(signed_polygon_area(R, Z));

            if areaAbs > bestArea
                bestArea = areaAbs;
                bestIdx = i;
            end
        end
    end

    if isempty(bestIdx)
        error('Could not select a usable flux contour.');
    end

    X = contours{bestIdx}.R(:);
    Z = contours{bestIdx}.Z(:);
end

% ======================================================================
function [X, Z] = clean_closed_curve(X, Z)
    X = X(:);
    Z = Z(:);

    % Remove repeated final point if present.
    if hypot(X(1)-X(end), Z(1)-Z(end)) < 1e-12
        X(end) = [];
        Z(end) = [];
    end

    % Remove consecutive duplicate points.
    d = hypot(diff([X; X(1)]), diff([Z; Z(1)]));
    keep = d(1:end-1) > 1e-12;

    X = X(keep);
    Z = Z(keep);

    if numel(X) < 20
        error('Contour has too few usable points after cleaning.');
    end
end

% ======================================================================
function A = signed_polygon_area(X, Z)
    X = X(:);
    Z = Z(:);

    X2 = [X(2:end); X(1)];
    Z2 = [Z(2:end); Z(1)];

    A = 0.5 * sum(X .* Z2 - X2 .* Z);
end

% ======================================================================
function [l, Xq, Zq, Ltot] = resample_closed_curve_by_arclength(X, Z, N)
    X = X(:);
    Z = Z(:);

    Xc = [X; X(1)];
    Zc = [Z; Z(1)];

    ds = hypot(diff(Xc), diff(Zc));
    s = [0; cumsum(ds)];
    Ltot = s(end);

    % Remove any zero-length segments before interpolation.
    [sUnique, ia] = unique(s, 'stable');
    Xc = Xc(ia);
    Zc = Zc(ia);

    l = linspace(0, Ltot, N+1).';
    l(end) = [];

    Xq = interp1(sUnique, Xc, l, 'pchip');
    Zq = interp1(sUnique, Zc, l, 'pchip');
end

% ======================================================================
function [ys, dydl, d2ydl2] = periodic_fourier_smooth_derivatives(y, Ltot, cutoff, p)
    y = y(:);
    N = numel(y);

    Y = fft(y);

    % Fourier mode numbers: 0,1,...,N/2,-N/2+1,...,-1
    if mod(N, 2) == 0
        m = [0:(N/2), -(N/2-1):-1].';
    else
        m = [0:((N-1)/2), -((N-1)/2):-1].';
    end

    cutoff = min(cutoff, floor(N/2)-1);

    % Smooth exponential low-pass filter.
    filt = exp(-(abs(m) / cutoff).^p);

    omega = 2*pi*m / Ltot;

    Ys = Y .* filt;

    ys      = real(ifft(Ys));
    dydl    = real(ifft(1i * omega .* Ys));
    d2ydl2  = real(ifft(-(omega.^2) .* Ys));
end

% ======================================================================
function plot_local_surface_result(surf)
    figure('Color', 'w', 'Name', 'Local flux surface construction');

    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % ------------------------------------------------------------
    % Raw vs smoothed surface
    % ------------------------------------------------------------
    nexttile;
    plot(surf.Xraw, surf.Zraw, '.', 'MarkerSize', 5); hold on;
    plot([surf.X; surf.X(1)], [surf.Z; surf.Z(1)], '-', 'LineWidth', 1.8);
    axis equal;
    grid on;
    xlabel('R or X');
    ylabel('Z');
    title(sprintf('Flux surface, \\psi_N = %.3f', surf.psin_target));
    legend('raw contour', 'smoothed periodic curve', 'Location', 'best');

    % ------------------------------------------------------------
    % Local coordinate mesh
    % ------------------------------------------------------------
    nexttile;
    plot([surf.X; surf.X(1)], [surf.Z; surf.Z(1)], 'k-', 'LineWidth', 1.5);
    hold on;

    % Plot several rho lines
    for j = 1:numel(surf.rho)
        plot(surf.Xrho(:, j), surf.Zrho(:, j), '-', 'LineWidth', 0.7);
    end

    % Plot a few normal lines
    idx = round(linspace(1, numel(surf.l), 24));
    rhoMin = min(surf.rho);
    rhoMax = max(surf.rho);

    for ii = idx
        xx = surf.X(ii) + surf.nX(ii) * [rhoMin rhoMax];
        zz = surf.Z(ii) + surf.nZ(ii) * [rhoMin rhoMax];
        plot(xx, zz, 'k-', 'LineWidth', 0.5);
    end

    axis equal;
    grid on;
    xlabel('R or X');
    ylabel('Z');
    title('Local coordinates: X(l,\rho), Z(l,\rho)');

    % ------------------------------------------------------------
    % u(l)
    % ------------------------------------------------------------
    nexttile;
    plot(surf.l, surf.u, 'LineWidth', 1.5);
    grid on;
    xlabel('l');
    ylabel('u(l)');
    title('Tangent angle');

    % ------------------------------------------------------------
    % Radius of curvature
    % ------------------------------------------------------------
    nexttile;
    plot(surf.l, surf.Rc_abs, 'LineWidth', 1.5);
    grid on;
    xlabel('l');
    ylabel('|R_c(l)|');
    title('Smoothed curvature radius');
end
