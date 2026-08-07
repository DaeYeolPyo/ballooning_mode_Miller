function ints = compute_contour_integrals(eq, eqfunc, surf)
    R = surf.R(:);
    Z = surf.Z(:);
    dl = surf.dl(:);

    val = eqfunc.eval(R, Z);

    gradPsi = val.gradPsi;
    F = val.F;

    if any(~isfinite(gradPsi)) || any(gradPsi <= 0)
        error('Invalid |grad psi| found on flux surface.');
    end

    % q(psi) from contour integral
    wq = F./(R.*gradPsi);
    q_int = sum(wq.*dl)/(2*pi);

    % GEQDSK q profile at this surface, for comparison
    q_geqdsk = eqfunc.q(surf.psiN);

    % Volume derivative
    wVprime = R./gradPsi;
    Vprime = 2*pi*sum(wVprime.*dl);

    % Enclosed toroidal volume from rotating the poloidal polygon
    [area_poloidal, intR_dA] = polygon_area_and_Rmoment(R, Z);
    V = 2*pi*abs(intR_dA);

    ints = struct();

    ints.psiN = surf.psiN;
    ints.psi = surf.psi;

    ints.q_int = q_int;
    ints.q_geqdsk = q_geqdsk;
    ints.q_relerr = (q_int - q_geqdsk)/max(1, abs(q_geqdsk));

    ints.V = V;
    ints.Vprime = Vprime;

    ints.area_poloidal = abs(area_poloidal);
    ints.mean_R = abs(intR_dA)/max(abs(area_poloidal), eps);
end

function [A, intR_dA] = polygon_area_and_Rmoment(R, Z)
    % For polygon vertices (R_i, Z_i),
    %   A = int dR dZ
    %   intR_dA = int R dR dZ
    
    R = R(:);
    Z = Z(:);

    Rp = [R(2:end); R(1)];
    Zp = [Z(2:end); Z(1)];

    cross = R.*Zp - Rp.*Z;

    A = 0.5*sum(cross);

    intR_dA = (1/6)*sum((R + Rp).*cross);
end