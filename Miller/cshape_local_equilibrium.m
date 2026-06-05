function base = cshape_local_equilibrium(eq, psin0, varargin)
%CSHAPE_LOCAL_EQUILIBRIUM Fit a GEQDSK surface to the C-shape local model.
%
%   base = cshape_local_equilibrium(eq, psin0)
%   base = cshape_local_equilibrium(filename, psin0, 'Name', value, ...)
%
% The fitted surface is
%   R = A + B sin(theta) + C cos(2 theta)
%   Z = G cos(theta) - H sin(2 theta)
%
% The output is a reusable one-surface geometry table for
% cshape_ballooning_coefficients.m.  The angular grid stored in base.theta
% is a PEST-like straight-field-line theta obtained from the local q
% integral, following the same remapping used by the Miller adapter.

    p = inputParser;
    addParameter(p, 'NTheta', 257, @(x)isnumeric(x)&&isscalar(x)&&x>=32);
    addParameter(p, 'NGeom', 1201, @(x)isnumeric(x)&&isscalar(x)&&x>=128);
    addParameter(p, 'Dpsin', 1e-4, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'UseBpFit', true, @(x)islogical(x)&&isscalar(x));
    addParameter(p, 'BN', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(p, 'aN', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
    addParameter(p, 'Q0', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)));
    addParameter(p, 'AlphaFactor', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x~=0));
    addParameter(p, 'CoeffOverride', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'InterpMethod', 'pchip', @(x)ischar(x)||isstring(x));
    parse(p, varargin{:});
    opt = p.Results;

    if ischar(eq) || isstring(eq)
        eq = read_geqdsk(char(eq));
    end

    ntheta = round(opt.NTheta);
    theta = linspace(0, 2*pi, ntheta+1).';
    theta(end) = [];

    if opt.UseBpFit
        coeffs = Cshape_metrics_bpfit(eq, psin0, round(opt.NGeom));
    else
        coeffs = Cshape_metrics(eq, psin0);
    end
    coeffs = apply_coeff_override(coeffs, opt.CoeffOverride);

    prof = cshape_profiles(eq, psin0, opt.InterpMethod);
    if isempty(opt.Q0)
        q0 = prof.q;
    else
        q0 = opt.Q0;
    end

    dpsi_dpsin = eq.sibry - eq.simag;
    if abs(dpsi_dpsin) < eps
        error('cshape_local_equilibrium:ZeroFluxSpan', ...
            'eq.sibry - eq.simag is zero.');
    end

    surf0 = cshape_surface_on_pest_theta(coeffs, theta, prof.F, ...
        dpsi_dpsin, round(opt.NGeom));

    coeffs_m = coeffs;
    coeffs_p = coeffs;
    coeffs_m(:,1) = coeffs(:,1) - opt.Dpsin .* coeffs(:,2);
    coeffs_p(:,1) = coeffs(:,1) + opt.Dpsin .* coeffs(:,2);

    surfm = cshape_surface_on_pest_theta(coeffs_m, theta, prof.F, ...
        dpsi_dpsin, round(opt.NGeom));
    surfp = cshape_surface_on_pest_theta(coeffs_p, theta, prof.F, ...
        dpsi_dpsin, round(opt.NGeom));

    R = surf0.R(:);
    Z = surf0.Z(:);
    R_psin = (surfp.R(:) - surfm.R(:)) ./ (2*opt.Dpsin);
    Z_psin = (surfp.Z(:) - surfm.Z(:)) ./ (2*opt.Dpsin);
    R_theta = periodic_derivative(theta, R);
    Z_theta = periodic_derivative(theta, Z);

    D = R_psin .* Z_theta - R_theta .* Z_psin;
    if any(abs(D) < 1e-12 * max(1, max(abs(D))))
        warning('cshape_local_equilibrium:SmallJacobian', ...
            'The fitted C-shape Jacobian has very small values.');
    end

    grad_psin  = [ Z_theta ./ D, zeros(ntheta,1), -R_theta ./ D ];
    grad_theta = [-Z_psin ./ D, zeros(ntheta,1),  R_psin ./ D ];
    grad_zeta  = [zeros(ntheta,1), 1./R, zeros(ntheta,1)];

    psi_R = dpsi_dpsin .* grad_psin(:,1);
    psi_Z = dpsi_dpsin .* grad_psin(:,3);
    BR = -psi_Z ./ R;
    BZ =  psi_R ./ R;
    Bphi = prof.F ./ R;
    Bvec = [BR, Bphi, BZ];
    Bmag = sqrt(sum(Bvec.^2, 2));

    if isempty(opt.BN)
        if isfield(eq, 'rmaxis') && isfinite(eq.rmaxis) && eq.rmaxis ~= 0
            BN = abs(prof.F ./ eq.rmaxis);
        elseif isfield(eq, 'bcentr') && isfinite(eq.bcentr) && eq.bcentr ~= 0
            BN = abs(eq.bcentr);
        else
            BN = max(Bmag);
        end
    else
        BN = opt.BN;
    end

    if isempty(opt.aN)
        aN = 0.5 .* (max(R) - min(R));
    else
        aN = opt.aN;
    end

    Bp2_psin = (surfp.Bp(:).^2 - surfm.Bp(:).^2) ./ (2*opt.Dpsin);

    alphaFactorDefault = default_alpha_factor(q0, eq, R, Bmag, ...
        dpsi_dpsin, grad_psin);
    if isempty(opt.AlphaFactor)
        alphaFactor = alphaFactorDefault;
    else
        alphaFactor = opt.AlphaFactor;
    end

    base = struct();
    base.model = 'cshape-local';
    base.psin = psin0;
    base.theta = theta;
    base.R = R;
    base.Z = Z;
    base.R_psin = R_psin;
    base.Z_psin = Z_psin;
    base.R_theta = R_theta;
    base.Z_theta = Z_theta;
    base.grad_psin = grad_psin;
    base.grad_theta = grad_theta;
    base.grad_zeta = grad_zeta;
    base.jacobian = D;
    base.BR = BR;
    base.Bphi = Bphi;
    base.BZ = BZ;
    base.B = Bmag;
    base.Bp = surf0.Bp(:);
    base.Bp2_psin = Bp2_psin;
    base.F = prof.F;
    base.q = q0;
    base.q_profile = prof.q;
    base.p = prof.p;
    base.pprime_geqdsk = prof.pprime_geqdsk;
    base.dpsi_dpsin = dpsi_dpsin;
    base.BN = BN;
    base.aN = aN;
    base.alphaFactor = alphaFactor;
    base.alphaFactorDefault = alphaFactorDefault;
    base.coeffs = coeffs;
    base.coeffNames = {'A','B','C','G','H'};
    base.surface = surf0;
    base.surface_minus = surfm;
    base.surface_plus = surfp;
    base.q_check = surf0.q_check;
    base.options = opt;
    base.eq = eq;
end

%==========================================================================
function coeffs = apply_coeff_override(coeffs, override)
    if isempty(override)
        return;
    end

    if isvector(override) && numel(override) == 5
        coeffs(:,1) = override(:);
        return;
    end

    if isequal(size(override), [5, 2])
        coeffs = override;
        return;
    end

    if isequal(size(override), [2, 5])
        coeffs = override.';
        return;
    end

    error('cshape_local_equilibrium:BadCoeffOverride', ...
        'CoeffOverride must be a 5-vector or a 5-by-2 coefficient table.');
end

%==========================================================================
function prof = cshape_profiles(eq, psin0, method)
    psin = linspace(0, 1, numel(eq.qpsi)).';
    method = char(method);

    prof = struct();
    prof.F = interp1(psin, eq.fpol(:), psin0, method, 'extrap');
    prof.q = interp1(psin, eq.qpsi(:), psin0, method, 'extrap');

    if isfield(eq, 'pres')
        prof.p = interp1(psin, eq.pres(:), psin0, method, 'extrap');
    else
        prof.p = NaN;
    end

    if isfield(eq, 'pprime')
        prof.pprime_geqdsk = interp1(psin, eq.pprime(:), psin0, method, 'extrap');
    else
        prof.pprime_geqdsk = NaN;
    end
end

%==========================================================================
function surf = cshape_surface_on_pest_theta(coeffs, theta_grid, F, dpsi_dpsin, ngeom)
    theta_raw = linspace(0, 2*pi, ngeom+1).';
    theta_raw(end) = [];

    [R, Z, R_theta, Z_theta, R_psin, Z_psin] = cshape_eval(coeffs, theta_raw);
    jac = R_psin .* Z_theta - R_theta .* Z_psin;
    htheta = hypot(R_theta, Z_theta);

    Bp = abs(dpsi_dpsin) .* htheta ./ max(abs(jac).*R, eps);
    pitch = abs(F) .* abs(jac) ./ max(abs(dpsi_dpsin).*R, eps);
    q_check = trapz_periodic(theta_raw, pitch) ./ (2*pi);

    theta_end = [theta_raw; 2*pi];
    pitch_end = [pitch; pitch(1)];
    theta_pest_end = [0; cumsum(0.5 .* ...
        (pitch_end(1:end-1) + pitch_end(2:end)) .* diff(theta_end))];

    if theta_pest_end(end) <= 0 || ~isfinite(theta_pest_end(end))
        error('cshape_local_equilibrium:BadPestMap', ...
            'Could not build a monotone PEST theta map.');
    end
    theta_pest = 2*pi .* theta_pest_end(1:end-1) ./ theta_pest_end(end);

    values = [R, Z, Bp, jac, htheta];
    values_q = interp_periodic_table(theta_pest, values, theta_grid(:));

    surf = struct();
    surf.theta = theta_grid(:);
    surf.R = values_q(:,1);
    surf.Z = values_q(:,2);
    surf.Bp = values_q(:,3);
    surf.jacobian = values_q(:,4);
    surf.htheta_geom = values_q(:,5);
    surf.theta_geom = theta_raw;
    surf.theta_pest = theta_pest;
    surf.q_check = q_check;
end

%==========================================================================
function [R, Z, R_theta, Z_theta, R_psin, Z_psin] = cshape_eval(coeffs, theta)
    c0 = coeffs(:,1);
    c1 = coeffs(:,2);

    A = c0(1); B = c0(2); C = c0(3); G = c0(4); H = c0(5);
    Ap = c1(1); Bp = c1(2); Cp = c1(3); Gp = c1(4); Hp = c1(5);

    theta = theta(:);
    R = A + B.*sin(theta) + C.*cos(2*theta);
    Z = G.*cos(theta) - H.*sin(2*theta);

    R_theta = B.*cos(theta) - 2*C.*sin(2*theta);
    Z_theta = -G.*sin(theta) - 2*H.*cos(2*theta);

    R_psin = Ap + Bp.*sin(theta) + Cp.*cos(2*theta);
    Z_psin = Gp.*cos(theta) - Hp.*sin(2*theta);
end

%==========================================================================
function alphaFactor = default_alpha_factor(q0, eq, R, Bmag, dpsi_dpsin, grad_psin)
    if isfield(eq, 'rmaxis') && isfinite(eq.rmaxis)
        R0 = eq.rmaxis;
    else
        R0 = mean(R, 'omitnan');
    end

    gradpsi = abs(dpsi_dpsin) .* vecnorm(grad_psin, 2, 2);
    gradpsiRef = mean(gradpsi, 'omitnan');
    B2ref = mean(Bmag.^2, 'omitnan');
    alphaFactor = -2 .* q0.^2 .* R0 .* gradpsiRef ./ max(B2ref, eps);
end

%==========================================================================
function df = periodic_derivative(x, f)
    x = x(:);
    f = f(:);
    h = mean(diff(x));
    df = (circshift(f,-1) - circshift(f,1)) ./ (2*h);
end

%==========================================================================
function val = trapz_periodic(theta, f)
    theta = theta(:);
    f = f(:);
    theta_end = [theta; 2*pi];
    f_end = [f; f(1)];
    val = sum(0.5 .* (f_end(1:end-1) + f_end(2:end)) .* diff(theta_end));
end

%==========================================================================
function vq = interp_periodic_table(theta, values, thetaq)
    theta = mod(theta(:), 2*pi);
    values = values(:,:);
    [theta, ord] = sort(theta);
    values = values(ord,:);
    [theta, keep] = unique(theta, 'stable');
    values = values(keep,:);

    theta_ext = [theta-2*pi; theta; theta+2*pi];
    values_ext = [values; values; values];
    thetaq = mod(thetaq(:), 2*pi);
    vq = interp1(theta_ext, values_ext, thetaq, 'pchip');
end
