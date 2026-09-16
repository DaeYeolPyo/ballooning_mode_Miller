function K = assemble_K_matrix(P3, quad, sym)
    arguments
        P3 (1,1) struct
        quad (1,1) struct
        sym logical = true
    end
    P1 = P3.P1points;
    T3 = P3.P1elements;

    T10 = P3.elements;

    Nt = size(T10, 1);
    Ndof = size(P3.points, 1);

    nEntry = Nt*100;
    I = zeros(nEntry, 1);
    J = zeros(nEntry, 1);
    V = zeros(nEntry, 1);

    counter = 0;

    for e = 1:Nt
        geomNodes = T3(e, :);
        Ve = P1(geomNodes, :);

        ids = T10(e, :);

        Ke = construct_K_matrix(Ve, quad, sym);

        [ii, jj] = ndgrid(ids, ids);

        idx = counter + (1:100);

        I(idx) = ii(:);
        J(idx) = jj(:);
        V(idx) = Ke(:);

        counter = counter + 100;
    end

    K = sparse(I, J, V, Ndof, Ndof);
end