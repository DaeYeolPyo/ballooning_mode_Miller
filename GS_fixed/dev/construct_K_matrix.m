function Ke = construct_K_matrix(Ve, quad, sym)
    arguments
        Ve (3,2) double
        quad (1,1) struct
        sym logical = false
    end

    %=============================================
    % Input check
    %=============================================
    if ~isequal(size(Ve), [3,2])
        error('[CONSTRUCT_K_MATRIX] Input vertices must be a 3 x 2 array.');
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

    a = J(1, 1);
    b = J(1, 2);
    c = J(2, 1);
    d = J(2, 2);

    invJT = (1/detJ)*[d, -c; -b, a];

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
    % Reference derivatives
    %
    % quad.dNdxi, quad.dNdeta : nq x 10
    %=============================================
    dNdxi  = quad.dNdxi;
    dNdeta = quad.dNdeta;

    % Physical gradients
    dNdR = invJT(1, 1)*dNdxi + invJT(1, 2)*dNdeta;
    dNdZ = invJT(2, 1)*dNdxi + invJT(2, 2)*dNdeta;

    % Quadrature coefficient
    alpha = quad.w.*detJ./Rq;

    %=============================================
    % Ke = dR' W dR + dZ' W dZ
    %=============================================
    Ke = dNdR.'*(alpha.*dNdR) + dNdZ.'*(alpha.*dNdZ);

    %=============================================
    % Symmetry cleanup
    %=============================================
    if sym
        Ke = 0.5*(Ke + Ke.');
    end
end