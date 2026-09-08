function f = assemble_RHS(P3, psi, quad, pprime, FFprime, dpsi)
    P1  = P3.P1points;
    T3  = P3.P1elements;
    T10 = P3.elements;

    Nt   = size(T10, 1);
    Ndof = size(P3.points, 1);

    psi = psi(:);

    if length(psi) ~= Ndof
        error('[ASSEMBLE_RHS] Global psi vector has incompatible size.');
    end

    f = zeros(Ndof, 1);

    %==================================================
    % Element assembly
    %==================================================
    for e = 1:Nt
        % Geometry vertices
        geomNodes = T3(e, :);
        Ve = P1(geomNodes, :);

        % P3 global DOFs
        ids = T10(e, :);

        % Local psi
        psi_e = psi(ids);

        % Local source
        fe = construct_RHS(Ve, psi_e, quad, pprime, FFprime, dpsi);

        % Assemble
        f(ids) = f(ids) + fe;
    end
end