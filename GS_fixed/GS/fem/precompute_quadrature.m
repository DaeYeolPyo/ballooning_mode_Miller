function quad = precompute_quadrature(p)
    arguments
        p (1,1) int8 = 5
    end

    [xi, eta, w] = Gauss_Dunavant_quadrature(p);

    nq = length(w);

    dNdxi  = zeros(nq, 10);
    dNdeta = zeros(nq, 10);
    N      = zeros(nq, 10);

    for q = 1:nq
        [Nq, dNref] = construct_P3_shape_functions(xi(q), eta(q));

        N(q, :)      = Nq;
        dNdxi(q, :)  = dNref(1, :);
        dNdeta(q, :) = dNref(2, :);
    end

    quad = struct();

    quad.xi  = xi(:);
    quad.eta = eta(:);
    quad.w   = w(:);

    quad.N      = N;
    quad.dNdxi  = dNdxi;
    quad.dNdeta = dNdeta;
end