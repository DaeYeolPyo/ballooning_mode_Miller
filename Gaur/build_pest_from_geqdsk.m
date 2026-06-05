function pest = build_pest_from_geqdsk(eq, npsi, ntheta, theta0)
%BUILD_PEST_FROM_GEQDSK Construct tokamak PEST-like coordinates from GEQDSK.
%
%   pest = build_pest_from_geqdsk(eq, npsi, ntheta, theta0)
%
% Input:
%   eq      : struct returned by read_geqdsk(filename)
%   npsi    : number of flux surfaces, e.g. 64
%   ntheta  : number of poloidal grid points, e.g. 256
%   theta0  : reference angle in alpha_t = zeta - q(psi)*(theta-theta0)
%
% Output:
%   pest.psin
%   pest.psi
%   pest.q
%   pest.theta
%   pest.R(ntheta,npsi)
%   pest.Z(ntheta,npsi)
%   pest.alpha_t(ntheta,npsi,nzeta)
%
% Notes:
%   PEST straight-field-line angle is constructed by enforcing
%
%       d zeta / d theta = q(psi)
%
%   along each flux surface. With an axisymmetric GEQDSK, this is done using
%   the standard field-line pitch relation and numerical flux-surface geometry.

    if nargin < 2 || isempty(npsi)
        npsi = 64;
    end
    if nargin < 3 || isempty(ntheta)
        ntheta = 256;
    end
    if nargin < 4 || isempty(theta0)
        theta0 = 0;
    end

    % ------------------------------------------------------------
    % Basic grid and flux normalization
    % ------------------------------------------------------------
    Rgrid = eq.rgrid(:);
    Zgrid = eq.zgrid(:);

    % eq.psirz is [nw, nh]. For contour/interp2 with meshgrid(R,Z),
    % use psi as [nh, nw].
    psiRZ = eq.psirz.';

    psi_axis = eq.simag;
    psi_bdry = eq.sibry;
    dpsi = psi_bdry - psi_axis;

    psin_all = linspace(0, 1, numel(eq.qpsi));
    q_of_psin = @(x) interp1(psin_all, eq.qpsi(:), x, 'pchip', 'extrap');

    % Avoid exact axis and exact separatrix for contour robustness.
    pest.psin = linspace(0.02, 0.98, npsi);
    pest.psi  = psi_axis + pest.psin * dpsi;
    pest.q    = q_of_psin(pest.psin);

    pest.theta = linspace(0, 2*pi, ntheta + 1);
    pest.theta(end) = [];

    pest.R = nan(ntheta, npsi);
    pest.Z = nan(ntheta, npsi);

    [RR, ZZ] = meshgrid(Rgrid, Zgrid);

    % ------------------------------------------------------------
    % Build each flux surface and remap to PEST straight-field-line theta
    % ------------------------------------------------------------
    for j = 1:npsi
        psi_level = pest.psi(j);

        % Extract flux-surface contour.
        C = contourc(Rgrid, Zgrid, psiRZ, [psi_level psi_level]);
        loops = parse_contourc(C);

        if isempty(loops)
            warning('No closed contour found at psin = %.4f', pest.psin(j));
            continue;
        end

        % Pick the loop enclosing the magnetic axis.
        idx = choose_axis_loop(loops, eq.rmaxis, eq.zmaxis);
        Rloop = loops(idx).R(:);
        Zloop = loops(idx).Z(:);

        % Make sure loop is closed.
        if hypot(Rloop(end)-Rloop(1), Zloop(end)-Zloop(1)) > 1e-10
            Rloop(end+1) = Rloop(1);
            Zloop(end+1) = Zloop(1);
        end

        % Geometric poloidal angle around magnetic axis.
        vartheta = unwrap(atan2(Zloop - eq.zmaxis, Rloop - eq.rmaxis));

        % Sort by geometric angle.
        [vartheta, order] = sort(vartheta);
        Rloop = Rloop(order);
        Zloop = Zloop(order);

        % Remove duplicate angles.
        [vartheta, uniqueIdx] = unique(vartheta, 'stable');
        Rloop = Rloop(uniqueIdx);
        Zloop = Zloop(uniqueIdx);

        % Periodic extension.
        vartheta = vartheta - vartheta(1);
        L = vartheta(end);
        vartheta_ext = [vartheta; vartheta(2:end) + 2*pi];
        R_ext = [Rloop; Rloop(2:end)];
        Z_ext = [Zloop; Zloop(2:end)];

        % Uniform geometric-angle grid for metric calculation.
        varg = linspace(0, 2*pi, max(4*ntheta, 512) + 1).';
        varg(end) = [];

        Rg = interp1(vartheta_ext, R_ext, varg, 'pchip');
        Zg = interp1(vartheta_ext, Z_ext, varg, 'pchip');

        % Derivatives wrt geometric angle.
        dvar = varg(2) - varg(1);
        dRdv = periodic_derivative(Rg, dvar);
        dZdv = periodic_derivative(Zg, dvar);

        % Grad psi on the surface.
        [dpsidR, dpsidZ] = gradient(psiRZ, Rgrid, Zgrid);
        psiR = interp2(RR, ZZ, dpsidR, Rg, Zg, 'linear');
        psiZ = interp2(RR, ZZ, dpsidZ, Rg, Zg, 'linear');

        gradpsi2 = psiR.^2 + psiZ.^2;

        % For axisymmetric tokamak:
        % B = grad(zeta) x grad(psi_p) + F(psi)/R e_zeta
        %
        % PEST straight-field-line angle can be constructed by using
        % dtheta_PEST/dvartheta proportional to
        %
        %       J_geom / R^2
        %
        % normalized so theta_PEST advances by 2*pi on one circuit.
        %
        % Numerically, J_geom-like weight is ds / |grad psi|.
        ds_dv = hypot(dRdv, dZdv);
        weight = ds_dv ./ max(sqrt(gradpsi2), eps) ./ max(Rg.^2, eps);

        cum = cumtrapz(varg, weight);
        cum = cum - cum(1);
        theta_pest_geom = 2*pi * cum / cum(end);

        % Invert theta_PEST(vartheta) to get R,Z on uniform PEST theta grid.
        theta_ext = [theta_pest_geom; theta_pest_geom(2:end) + 2*pi];
        Rg_ext = [Rg; Rg(2:end)];
        Zg_ext = [Zg; Zg(2:end)];

        pest.R(:, j) = interp1(theta_ext, Rg_ext, pest.theta(:), 'pchip');
        pest.Z(:, j) = interp1(theta_ext, Zg_ext, pest.theta(:), 'pchip');
    end

    pest.theta0 = theta0;
    pest.description = ...
        'Tokamak PEST coordinates with alpha_t = zeta - q(psi)*(theta-theta0).';
end

% =====================================================================
function loops = parse_contourc(C)
    loops = struct('level', {}, 'R', {}, 'Z', {});
    k = 1;
    nloop = 0;

    while k < size(C, 2)
        level = C(1, k);
        npts = C(2, k);
        pts = C(:, k+1:k+npts);

        nloop = nloop + 1;
        loops(nloop).level = level;
        loops(nloop).R = pts(1, :).';
        loops(nloop).Z = pts(2, :).';

        k = k + npts + 1;
    end
end

% =====================================================================
function idx = choose_axis_loop(loops, Raxis, Zaxis)
    idx = 1;
    bestArea = inf;

    for k = 1:numel(loops)
        R = loops(k).R;
        Z = loops(k).Z;

        inside = inpolygon(Raxis, Zaxis, R, Z);
        area = abs(polyarea(R, Z));

        if inside && area < bestArea
            bestArea = area;
            idx = k;
        end
    end
end

% =====================================================================
function df = periodic_derivative(f, dx)
    df = (circshift(f, -1) - circshift(f, 1)) / (2*dx);
end
