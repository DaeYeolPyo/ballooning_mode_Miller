function bln = compute_bishop_ch3_ballooning_terms(surf, opts)
%COMPUTE_BISHOP_CH3_BALLOONING_TERMS
%
% Compute Bishop Chapter 3 ballooning-equation coefficients on one
% selected GEQDSK flux surface.
%
% Input:
%   surf : output of build_local_flux_surface_from_geqdsk(...)
%
% Required fields in surf:
%   surf.eq
%   surf.l
%   surf.X            major radius X(l)
%   surf.Z
%   surf.u            tangent angle, dX/dl = cos u, dZ/dl = sin u
%   surf.Rc_signed or surf.Rc_abs
%   surf.psin_target
%
% Output:
%   bln.Bp            B_p^(0)(l)
%   bln.B2            B^2 = Bp^2 + I^2/X^2
%   bln.gradS2        Eq. (29), |grad S|^2
%   bln.curvTerm      Eq. (31), 2 p' (B x gradS).kappa / B^2
%   bln.A             coefficient inside derivative in Eq. (32):
%                     A = (|gradS|^2 / B^2) Bp
%   bln.C             curvature coefficient in Eq. (32):
%                     C = Eq. (31)
%   bln.L0            if F is supplied later:
%                     Bp d/dl [ A dF/dl ] + C F = 0
%
% Notes:
%   - This follows the notation of Bishop Chapter 3.
%   - In Bishop's notation x0 h0 = X(l).
%   - I is the toroidal field function, I = R B_phi = fpol.
%   - Iprime is dI/dpsi, computed from ffprim/I.
%   - pprime is dp/dpsi.
%
% Usage:
%   bln = compute_bishop_ch3_ballooning_terms(surf);
%
%   figure; plot(bln.l, bln.Bp); xlabel('\ell'); ylabel('B_p');
%   figure; plot(bln.l, bln.gradS2); xlabel('\ell'); ylabel('|grad S|^2');
%   figure; plot(bln.l, bln.curvTerm); xlabel('\ell'); ylabel('curvature term');

    if nargin < 2
        opts = struct();
    end
    opts = set_default_ballooning_opts(opts);

    eq = surf.eq;

    l = surf.l(:);
    X = surf.X(:);
    Z = surf.Z(:);
    u = surf.u(:);

    N = numel(l);
    L = surf.L;

    if isfield(surf, 'Rc_signed')
        Rcurv = surf.Rc_signed(:);
    else
        Rcurv = surf.Rc_abs(:);
    end

    % Avoid division by exactly zero curvature radius.
    Rcurv(abs(Rcurv) < opts.minCurvatureRadius) = ...
        sign(Rcurv(abs(Rcurv) < opts.minCurvatureRadius)) .* opts.minCurvatureRadius;

    % ------------------------------------------------------------------
    % 0. Choose x0 and h0
    % ------------------------------------------------------------------
    % Bishop uses x = x0 h0. Since X(l) is the actual cylindrical
    % major radius on the flux surface, we can set
    %
    %   h0(l) = X(l) / x0.
    %
    % The origin l=0 is arbitrary, so x0 is the major radius at the
    % reference index. By default we use index 1.
    if isempty(opts.refIndex)
        i0 = 1;
    else
        i0 = opts.refIndex;
    end

    x0 = X(i0);
    h0 = X ./ x0;

    % ------------------------------------------------------------------
    % 1. Compute B_p(l) from GEQDSK psi(R,Z)
    % ------------------------------------------------------------------
    [Bp, dpsi_dR_s, dpsi_dZ_s] = compute_Bp_on_surface(eq, X, Z, opts);

    % Optional smoothing of Bp along l.
    if opts.smoothBp
        Bp = periodic_smooth_real(Bp, opts.BpSmoothCutoff, opts.filterPower);
    end

    Bp = max(Bp, opts.minBp);

    % ------------------------------------------------------------------
    % 2. Surface quantities I, I', p'
    % ------------------------------------------------------------------
    psin = eq.psin(:);
    fpol = eq.fpol(:);
    pprime_prof = eq.pprime(:);
    ffprim_prof = eq.ffprim(:);

    psin0 = surf.psin_target;

    I0 = interp1(psin, fpol, psin0, opts.profileInterp, 'extrap');

    ffprim0 = interp1(psin, ffprim_prof, psin0, opts.profileInterp, 'extrap');
    pprime0 = interp1(psin, pprime_prof, psin0, opts.profileInterp, 'extrap');

    % GEQDSK ffprim = F F', with F = fpol = I.
    Iprime0 = ffprim0 / I0;

    % Allow manual override if user wants to scan eigenvalue pprime.
    if ~isempty(opts.I0)
        I0 = opts.I0;
    end
    if ~isempty(opts.Iprime0)
        Iprime0 = opts.Iprime0;
    end
    if ~isempty(opts.pprime0)
        pprime0 = opts.pprime0;
    end

    % ------------------------------------------------------------------
    % 3. Useful magnetic quantities
    % ------------------------------------------------------------------
    Xgeom = x0 .* h0;  % should equal X

    Btor = I0 ./ Xgeom;
    B2 = Bp.^2 + Btor.^2;
    B4 = B2.^2;

    dB2dl = periodic_derivative_real(B2, L);
    if opts.smoothDerivatives
        dB2dl = periodic_smooth_real(dB2dl, opts.derivSmoothCutoff, opts.filterPower);
    end

    % ------------------------------------------------------------------
    % 4. Integral appearing in Eq. (29) and Eq. (31)
    % ------------------------------------------------------------------
    %
    % J(l) = integral_{l0}^{l} dl/(h0^2 Bp) * Q(l)
    %
    % where
    %
    % Q = I'/I * (1 + I^2/(x0^2 h0^2 Bp^2))
    %     + p'/Bp^2
    %     + 2/R * (1 + R sin(u)/(x0 h0))
    %     - 1/(x0 h0 Bp)
    %
    % This is the bracketed integral in Eq. (29).
    %
    Q29 = ...
        (Iprime0 / I0) .* (1 + I0^2 ./ (x0^2 .* h0.^2 .* Bp.^2)) ...
        + pprime0 ./ Bp.^2 ...
        + (2 ./ Rcurv) .* (1 + Rcurv .* sin(u) ./ (x0 .* h0)) ...
        - 1 ./ (x0 .* h0 .* Bp);

    integrandJ = Q29 ./ (h0.^2 .* Bp);

    J = cumulative_trapezoid_periodic_path(l, integrandJ);

    % ------------------------------------------------------------------
    % 5. Eq. (29): |grad S|^2
    % ------------------------------------------------------------------
    gradS2_first = ...
        1 ./ (x0^2 .* h0.^2) .* ...
        (1 + I0^2 ./ (x0^2 .* h0.^2 .* Bp.^2));

    gradS2_shear = ...
        (h0.^2 .* Bp.^2 .* I0^2 ./ x0^2) .* J.^2;

    gradS2 = gradS2_first + gradS2_shear;

    % ------------------------------------------------------------------
    % 6. Eq. (31): curvature drive term
    % ------------------------------------------------------------------
    %
    % Eq. (31):
    %
    % 2 p' (B x gradS).kappa / B^2
    %   = 2*p'/(x0 h0) * {
    %       2 I^2 sin u / (x0^3 h0^3 B^2 Bp)
    %       - 2 Bp / (R B^2)
    %       - I^2 h0 Bp/(x0 B^4) * dB^2/dl * J
    %     }
    %
    % where J is the same integral defined above.
    %
    curv_piece_toroidicity = ...
        2 .* I0^2 .* sin(u) ./ ...
        (x0^3 .* h0.^3 .* B2 .* Bp);

    curv_piece_normalcurv = ...
        -2 .* Bp ./ (Rcurv .* B2);

    curv_piece_shear = ...
        - I0^2 .* h0 .* Bp ./ (x0 .* B4) .* dB2dl .* J;

    curvBracket = curv_piece_toroidicity + ...
                  curv_piece_normalcurv + ...
                  curv_piece_shear;

    curvTerm = 2 .* pprime0 ./ (x0 .* h0) .* curvBracket;

    % ------------------------------------------------------------------
    % 7. Eq. (32): ballooning equation coefficients
    % ------------------------------------------------------------------
    %
    % Bishop Eq. (32):
    %
    %   Bp d/dl { (|gradS|^2/B^2) Bp dF/dl }
    %   + curvTerm F = 0
    %
    % Define
    %
    %   A(l) = (|gradS|^2/B^2) Bp
    %   C(l) = curvTerm
    %
    A = gradS2 ./ B2 .* Bp;
    C = curvTerm;

    % ------------------------------------------------------------------
    % 8. Pack output
    % ------------------------------------------------------------------
    bln = struct();

    bln.l = l;
    bln.L = L;

    bln.X = X;
    bln.Z = Z;
    bln.x0 = x0;
    bln.h0 = h0;
    bln.u = u;
    bln.Rcurv = Rcurv;

    bln.Bp = Bp;
    bln.Btor = Btor;
    bln.B2 = B2;
    bln.B4 = B4;
    bln.dB2dl = dB2dl;

    bln.I0 = I0;
    bln.Iprime0 = Iprime0;
    bln.pprime0 = pprime0;

    bln.dpsi_dR = dpsi_dR_s;
    bln.dpsi_dZ = dpsi_dZ_s;

    bln.Q29 = Q29;
    bln.integrandJ = integrandJ;
    bln.J = J;

    bln.gradS2_first = gradS2_first;
    bln.gradS2_shear = gradS2_shear;
    bln.gradS2 = gradS2;

    bln.curv_piece_toroidicity = curv_piece_toroidicity;
    bln.curv_piece_normalcurv = curv_piece_normalcurv;
    bln.curv_piece_shear = curv_piece_shear;
    bln.curvBracket = curvBracket;
    bln.curvTerm = curvTerm;

    bln.A = A;
    bln.C = C;

    bln.opts = opts;

    if opts.plot
        plot_bishop_ch3_terms(bln);
    end
end

% ======================================================================
function opts = set_default_ballooning_opts(opts)

    if ~isfield(opts, 'refIndex')
        opts.refIndex = [];
    end

    if ~isfield(opts, 'profileInterp')
        opts.profileInterp = 'pchip';
    end

    if ~isfield(opts, 'I0')
        opts.I0 = [];
    end

    if ~isfield(opts, 'Iprime0')
        opts.Iprime0 = [];
    end

    if ~isfield(opts, 'pprime0')
        opts.pprime0 = [];
    end

    if ~isfield(opts, 'smoothBp')
        opts.smoothBp = true;
    end

    if ~isfield(opts, 'smoothDerivatives')
        opts.smoothDerivatives = true;
    end

    if ~isfield(opts, 'BpSmoothCutoff')
        opts.BpSmoothCutoff = 80;
    end

    if ~isfield(opts, 'derivSmoothCutoff')
        opts.derivSmoothCutoff = 80;
    end

    if ~isfield(opts, 'filterPower')
        opts.filterPower = 8;
    end

    if ~isfield(opts, 'minBp')
        opts.minBp = 1e-10;
    end

    if ~isfield(opts, 'minCurvatureRadius')
        opts.minCurvatureRadius = 1e-8;
    end

    if ~isfield(opts, 'plot')
        opts.plot = true;
    end
end

% ======================================================================
function [Bp, dpsi_dR_s, dpsi_dZ_s] = compute_Bp_on_surface(eq, X, Z, opts)

    Rgrid = eq.rgrid(:).';
    Zgrid = eq.zgrid(:);

    psiRZ = eq.psirz.';  % size: [nh, nw], rows Z, columns R

    dR = Rgrid(2) - Rgrid(1);
    dZ = Zgrid(2) - Zgrid(1);

    [dpsi_dZ_grid, dpsi_dR_grid] = gradient(psiRZ, dZ, dR);

    % Use griddedInterpolant instead of interp2 for cleaner extrap handling.
    FR = griddedInterpolant({Zgrid, Rgrid}, dpsi_dR_grid, 'spline', 'nearest');
    FZ = griddedInterpolant({Zgrid, Rgrid}, dpsi_dZ_grid, 'spline', 'nearest');

    dpsi_dR_s = FR(Z, X);
    dpsi_dZ_s = FZ(Z, X);

    % Axisymmetric poloidal field:
    %   B_R = -1/R dpsi/dZ
    %   B_Z =  1/R dpsi/dR
    %
    % Therefore:
    %   Bp = |grad psi| / R.
    Bp = sqrt(dpsi_dR_s.^2 + dpsi_dZ_s.^2) ./ X;

    if opts.smoothBp
        Bp = periodic_smooth_real(Bp, opts.BpSmoothCutoff, opts.filterPower);
    end
end

% ======================================================================
function J = cumulative_trapezoid_periodic_path(l, f)
    l = l(:);
    f = f(:);

    N = numel(l);
    J = zeros(N, 1);

    dl = diff(l);

    J(2:end) = cumsum(0.5 .* (f(1:end-1) + f(2:end)) .* dl);
end

% ======================================================================
function dydl = periodic_derivative_real(y, L)
    y = y(:);
    N = numel(y);

    Y = fft(y);

    if mod(N, 2) == 0
        m = [0:(N/2), -(N/2-1):-1].';
    else
        m = [0:((N-1)/2), -((N-1)/2):-1].';
    end

    omega = 2*pi*m/L;

    dydl = real(ifft(1i .* omega .* Y));
end

% ======================================================================
function ys = periodic_smooth_real(y, cutoff, p)
    y = y(:);
    N = numel(y);

    Y = fft(y);

    if mod(N, 2) == 0
        m = [0:(N/2), -(N/2-1):-1].';
    else
        m = [0:((N-1)/2), -((N-1)/2):-1].';
    end

    cutoff = min(cutoff, floor(N/2)-1);
    filt = exp(-(abs(m) ./ cutoff).^p);

    ys = real(ifft(Y .* filt));
end

% ======================================================================
function plot_bishop_ch3_terms(bln)

    figure('Color', 'w', 'Name', 'Bishop Chapter 3 ballooning terms');

    tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(bln.l, bln.Bp, 'LineWidth', 1.5);
    grid on;
    xlabel('\ell');
    ylabel('B_p^{(0)}');
    title('Poloidal field');

    nexttile;
    plot(bln.l, bln.B2, 'LineWidth', 1.5);
    grid on;
    xlabel('\ell');
    ylabel('B^2');
    title('Total magnetic field squared');

    nexttile;
    plot(bln.l, bln.gradS2, 'LineWidth', 1.5);
    hold on;
    plot(bln.l, bln.gradS2_first, '--', 'LineWidth', 1.0);
    plot(bln.l, bln.gradS2_shear, ':', 'LineWidth', 1.2);
    grid on;
    xlabel('\ell');
    ylabel('| \nabla S |^2');
    title('Eq. (29)');
    legend('total', 'metric part', 'integral/shear part', 'Location', 'best');

    nexttile;
    plot(bln.l, bln.J, 'LineWidth', 1.5);
    grid on;
    xlabel('\ell');
    ylabel('J(\ell)');
    title('Integral in Eq. (29), Eq. (31)');

    nexttile;
    plot(bln.l, bln.curvTerm, 'LineWidth', 1.5);
    grid on;
    xlabel('\ell');
    ylabel('curvature term');
    title('Eq. (31)');

    nexttile;
    plot(bln.l, bln.A, 'LineWidth', 1.5);
    hold on;
    plot(bln.l, bln.C, 'LineWidth', 1.5);
    grid on;
    xlabel('\ell');
    ylabel('coefficient');
    title('Eq. (32) coefficients');
    legend('A=(|\nabla S|^2/B^2)B_p', 'C=curvature term', ...
           'Location', 'best');
end
