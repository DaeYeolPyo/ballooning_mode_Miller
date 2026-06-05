function bal = cshape_ballooning_coefficients(base, varargin)
%CSHAPE_BALLOONING_COEFFICIENTS Build Gaur-style coefficients for C-shape.
%
%   bal = cshape_ballooning_coefficients(base, 'S_hat', s, 'Alpha', alpha)
%
% Input
%   base : output of cshape_local_equilibrium.m
%
% The output struct is compatible with ../Gaur/solve_ballooning_eigenvalue.m.
% The default shear definition is the usual flux-coordinate estimate
%
%   s_hat = 2*psi_N/q * dq/dpsi_N.
%
% Use 'ShearDefinition','minor' to instead use the Miller adapter convention
% dq/dpsi_N = q*s_hat.

    p = inputParser;
    addParameter(p, 'S_hat', 0.0, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Alpha', 0.0, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'Theta0', 0.0, @(x)isnumeric(x)&&isscalar(x));
    addParameter(p, 'ShearDefinition', 'flux', @(x)ischar(x)||isstring(x));
    addParameter(p, 'AlphaFactor', [], @(x) isempty(x) || (isnumeric(x)&&isscalar(x)&&x~=0));
    addParameter(p, 'IncludeFPrime', true, @(x)islogical(x)&&isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    theta = base.theta(:);
    R = base.R(:);

    q0 = base.q;
    F0 = base.F;
    BN = base.BN;
    aN = base.aN;
    dpsi_dpsin = base.dpsi_dpsin;

    dq_dpsin = shear_to_dq_dpsin(opt.S_hat, q0, base.psin, opt.ShearDefinition);

    if isempty(opt.AlphaFactor)
        alphaFactor = base.alphaFactor;
    else
        alphaFactor = opt.AlphaFactor;
    end

    % alpha = alphaFactor * (mu0 dp/dpsi_GEQDSK)
    pprimeEquation = opt.Alpha ./ alphaFactor;
    dpbar_dpsin = pprimeEquation .* dpsi_dpsin ./ (BN.^2);

    if opt.IncludeFPrime && abs(q0) > eps
        F_psin = F0 .* dq_dpsin ./ q0;
    else
        F_psin = 0;
    end

    grad_psin = base.grad_psin;
    grad_theta = base.grad_theta;
    grad_zeta = base.grad_zeta;

    psi_R = dpsi_dpsin .* grad_psin(:,1);
    psi_Z = dpsi_dpsin .* grad_psin(:,3);
    BR = -psi_Z ./ R;
    BZ =  psi_R ./ R;
    Bphi = F0 ./ R;
    Bvec = [BR, Bphi, BZ];
    Bmag = sqrt(sum(Bvec.^2, 2));
    bvec = Bvec ./ Bmag;
    BoverBN = Bmag ./ BN;

    B2 = Bmag.^2;
    B2_theta = periodic_derivative(theta, B2);
    Bphi2_psin = 2.*F0.*F_psin./(R.^2) ...
        - 2.*F0.^2.*base.R_psin(:)./(R.^3);
    B2_psin = base.Bp2_psin(:) + Bphi2_psin;
    grad_B2 = B2_psin .* grad_psin + B2_theta .* grad_theta;

    grad_H = dpbar_dpsin .* grad_psin + (0.5 ./ BN.^2) .* grad_B2;

    gradN_psin = aN .* grad_psin;
    gradN_theta = aN .* grad_theta;
    gradN_zeta = aN .* grad_zeta;
    gradN_H = aN .* grad_H;

    theta_shift = theta - opt.Theta0;
    gradN_alpha = gradN_zeta ...
        - q0 .* gradN_theta ...
        - dq_dpsin .* theta_shift .* gradN_psin;

    b_dot_gradN_theta = sum(bvec .* gradN_theta, 2);
    if mean(b_dot_gradN_theta, 'omitnan') < 0
        warning('cshape_ballooning_coefficients:NegativeBdotGradTheta', ...
            ['b dot grad_N theta is negative on average. The fitted theta ', ...
             'orientation may be opposite to the field-line convention.']);
    end

    gradN_alpha_sq = sum(gradN_alpha.^2, 2);
    cross_term = sum(cross(bvec, gradN_H, 2) .* gradN_alpha, 2);

    gcoef = b_dot_gradN_theta .* gradN_alpha_sq ./ BoverBN;
    fcoef = gradN_alpha_sq ./ (b_dot_gradN_theta .* BoverBN.^3);
    ccoef = 2 ./ b_dot_gradN_theta ./ BoverBN.^4 ...
        .* dpbar_dpsin .* cross_term;

    bal = struct();
    bal.model = 'cshape-local';
    bal.psin = base.psin;
    bal.theta = theta;
    bal.R = R;
    bal.Z = base.Z(:);
    bal.g = gcoef;
    bal.c = ccoef;
    bal.f = fcoef;
    bal.q = q0;
    bal.dq_dpsin = dq_dpsin;
    bal.s_hat = opt.S_hat;
    bal.alpha = opt.Alpha;
    bal.alphaFactor = alphaFactor;
    bal.pprimeEquation = pprimeEquation;
    bal.dpbar_dpsin = dpbar_dpsin;
    bal.F = F0;
    bal.F_psin = F_psin;
    bal.BN = BN;
    bal.aN = aN;
    bal.theta0 = opt.Theta0;
    bal.shearDefinition = char(lower(string(opt.ShearDefinition)));

    bal.BR = BR;
    bal.Bphi = Bphi;
    bal.BZ = BZ;
    bal.bvec = bvec;
    bal.B = Bmag;
    bal.BoverBN = BoverBN;
    bal.b_dot_gradN_theta = b_dot_gradN_theta;
    bal.gradN_alpha_sq = gradN_alpha_sq;
    bal.gradN_alpha = gradN_alpha;
    bal.gradN_psin = gradN_psin;
    bal.gradN_theta = gradN_theta;
    bal.gradN_zeta = gradN_zeta;
    bal.gradN_H = gradN_H;
    bal.cross_term = cross_term;
    bal.B2_psin = B2_psin;
    bal.B2_theta = B2_theta;
    bal.cshape = base;
    bal.options = opt;
end

%==========================================================================
function dq = shear_to_dq_dpsin(sHat, q0, psin0, definition)
    definition = lower(string(definition));
    switch definition
        case "flux"
            dq = sHat .* q0 ./ max(2.*psin0, eps);
        case "minor"
            dq = sHat .* q0;
        otherwise
            error('cshape_ballooning_coefficients:BadShearDefinition', ...
                'ShearDefinition must be ''flux'' or ''minor''.');
    end
end

%==========================================================================
function df = periodic_derivative(x, f)
    x = x(:);
    f = f(:);
    h = mean(diff(x));
    df = (circshift(f,-1) - circshift(f,1)) ./ (2*h);
end
