function fe = construct_RHS(Ve, psi_e, quad, pprime, FFprime, dpsi)
% [CONSTRUCT_RHS]
%
% Construct the 10 x 1 local source vector for the P3
% Grad-Shafranov finite element.
%
% Weak form:
%
%   fe(i) = integral_Omega_e N_i * S(R,psi) dR dZ
%
% where
%
%   S(R,psi)
%       = mu0 * R * p'(psi)
%       + FF'(psi) / R
%
%
% INPUT
%
%   Ve : 3 x 2
%       Physical coordinates of the affine triangle vertices:
%
%           [R1 Z1
%            R2 Z2
%            R3 Z3]
%
%   psi_e : 10 x 1
%       P3 nodal values of psi in this element.
%
%   quad : structure containing
%
%       quad.xi      : nq x 1
%       quad.eta     : nq x 1
%       quad.w       : nq x 1
%       quad.N       : nq x 10
%
%   pprime
%       array containing dp/dpsi_N
%
%   FFprime
%       array containing F*dF/dpsi_N
%
%   dpsi
%       psi_boundary - psi_axis
%
% OUTPUT
%
%   fe : 10 x 1 local source vector
%
    mu0 = 4*pi*1.e-7;

    %=====================================================
    % Input checks
    %=====================================================
    if ~isequal(size(Ve), [3,2])
        error('[CONSTRUCT_RHS] Ve must be a 3 x 2 array.');
    end

    psi_e = psi_e(:);

    if length(psi_e) ~= 10
        error('[CONSTRUCT_RHS] psi_e must contain 10 P3 nodal value.');
    end

    %=============================================
    % Jacobian
    %=============================================
    R1 = Ve(1,1);  Z1 = Ve(1,2);
    R2 = Ve(2,1);  Z2 = Ve(2,2);
    R3 = Ve(3,1);  Z3 = Ve(3,2);

    J = [R2-R1, R3-R1; Z2-Z1, Z3-Z1];

    detJ = det(J);

    if detJ <= 0
        error('[CONSTRUCT_K_MATRIX] Invalid element: det(J) <= 0.');
    end

    %=============================================
    % R evaluation
    %=============================================
    L1 = 1 - quad.xi - quad.eta;
    L2 = quad.xi;
    L3 = quad.eta;

    Rq = L1*R1 + L2*R2 + L3*R3;

    if any(Rq <= 0)
        error('[CONSTRUCT_K_MATRIX] Encountered R <= 0.')
    end

    %=============================================
    % psi at quadrature points
    %
    % quad.N : nq x 10
    % psi_e  : 10 x 1
    %
    % psi_q  : nq x 1
    %=============================================
    psi_q = quad.N*psi_e;

    %=============================================
    % Evaluate profiles
    %=============================================
    nw = length(pprime);
    psiN = linspace(0, 1, nw);

    pprime_func  = griddedInterpolant(psiN, pprime/dpsi);
    FFprime_func = griddedInterpolant(psiN, FFprime/dpsi);

    % For psi_boundary = 0,
    % psiN = (psi - psi_axis)/(-psi_axis) where -psi_axis = dpsi
    psiN_q = (psi_q + dpsi)/dpsi;

    pprime_q  = pprime_func(psiN_q);
    FFprime_q = FFprime_func(psiN_q);

    %=============================================
    % RHS
    %
    % mu0 R p' + FF'/R
    %=============================================
    Sq = mu0.*Rq.*pprime_q + FFprime_q./Rq;

    %=============================================
    % Quadrature weights
    %=============================================
    alpha = quad.w.*detJ;

    %=============================================
    % Local load vector
    %
    % fe_i = sum_q N_i(q) * S_q * w_q * detJ
    %
    % Matrix form:
    %
    % fe = N' * (alpha .* S)
    %=============================================
    fe = quad.N.'*(alpha.*Sq);
end