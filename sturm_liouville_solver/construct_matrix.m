function [Gmat, Cmat, Fmat, thetaDof] = construct_matrix(thetaGrid, g, c, f)
    % Convert all inputs to column vectors
    thetaGrid = thetaGrid(:);
    g = g(:);
    c = c(:);
    f = f(:);

    ngrid = numel(thetaGrid);

    if ~isvector(thetaGrid) || ...
            numel(g) ~= ngrid || ...
            numel(c) ~= ngrid || ...
            numel(f) ~= ngrid
        error(['thetaGrid, g, c, and f must be vectors ', ...
            'with the same number of entries.']);
    end

    if any(diff(thetaGrid) <= 0)
        error('thetaGrid must be strictly increasing.');
    end

    if any(~isfinite(thetaGrid)) || ...
            any(~isfinite(g)) || ...
            any(~isfinite(c)) || ...
            any(~isfinite(f))
        error('thetaGrid, g, c, and f must contain finite values.');
    end

    if ngrid < 4
        error('At least four grid points are required for interpolation.');
    end

    % thetaGrid contains element boundaries
    nelm = ngrid - 1;

    % P2 Lagrange basis: endpoint DOFs + one midpoint DOF per element
    ndof = 2*nelm + 1;

    % Coordinates corresponding to the global P2 DOFs
    thetaDof = zeros(ndof, 1);
    thetaDof(1:2:end) = thetaGrid;
    thetaDof(2:2:end) = 0.5*(thetaGrid(1:end-1) + thetaGrid(2:end));

    % Global matrices
    % Each P2 element contributes a 3-by-3 block.
    Gmat = spalloc(ndof, ndof, 9*nelm);
    Cmat = spalloc(ndof, ndof, 9*nelm);
    Fmat = spalloc(ndof, ndof, 9*nelm);

    % Spline interpolation
    gfunc = griddedInterpolant(thetaGrid, g, 'pchip');
    cfunc = griddedInterpolant(thetaGrid, c, 'pchip');
    ffunc = griddedInterpolant(thetaGrid, f, 'pchip');

    % Using nq-point Gaussian-Lagrange quadrature
    % for integration within elements.
    nq = 4;
    [xi_q, w_q] = GauLeg_points(nq);

    xi_q = xi_q(:);
    w_q = w_q(:);

    % Pre-compute reference element P2 basis functions
    Nq = zeros(nq, 3);
    dNq_dxi = zeros(nq, 3);

    for q = 1:nq
        xi = xi_q(q);

        % P2 shape functions
        % N1 = (1/2)xi(xi - 1)
        % N2 = 1 - xi^2
        % N3 = (1/2)xi(xi + 1)
        Nq(q, :) = [0.5*xi*(xi - 1), 1 - xi^2, 0.5*xi*(xi + 1)];
        dNq_dxi(q, :) = [xi - 0.5, -2*xi, xi + 0.5];
    end

    % Element loop
    for e = 1:nelm
        % P2 local-to-global connectivity
        I = [2*e - 1, 2*e, 2*e + 1];

        thL = thetaGrid(e);
        thR = thetaGrid(e+1);

        he = thR - thL;
        Je = 0.5*he;

        Ge = zeros(3, 3);
        Ce = zeros(3, 3);
        Fe = zeros(3, 3);

        % Gaussian quadrature
        for q = 1:nq
            % Gaussian abscissa and weights
            xi = xi_q(q);
            w = w_q(q);

            Nshape = Nq(q, :);
            dN_dxi = dNq_dxi(q, :);

            % Physical quadrature point
            thq = 0.5*(thL + thR) + 0.5*he*xi;

            % Physical-coordinate derivatives
            dN_dth = dN_dxi/Je;

            % Coefficients evaluated in physical coordinates
            gq = gfunc(thq);
            cq = cfunc(thq);
            fq = ffunc(thq);

            % Element matrices
            Ge = Ge + Je*w*gq*(dN_dth.'*dN_dth);
            Ce = Ce + Je*w*cq*(Nshape.'*Nshape);
            Fe = Fe + Je*w*fq*(Nshape.'*Nshape);
        end

        % Global assembly
        Gmat(I, I) = Gmat(I, I) + Ge;
        Cmat(I, I) = Cmat(I, I) + Ce;
        Fmat(I, I) = Fmat(I, I) + Fe;
    end

    % Remove tiny asymmetries caused by roundoff
%     Gmat = 0.5*(Gmat + Gmat.');
%     Cmat = 0.5*(Cmat + Cmat.');
%     Fmat = 0.5*(Fmat + Fmat.');
end