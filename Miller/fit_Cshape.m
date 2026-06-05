function coeffs = fit_Cshape(surf)
    R_Z0 = find_R_at_Z0_from_surface(surf);

    Rmid = mean(R_Z0);
    dR = abs(R_Z0(2) - R_Z0(1));
    Rmax = max(surf.R);
    [Zmax, thZInd] = max(surf.Z);
    Rtop = surf.R(thZInd);

    D = Rmax - Rmid;

    B = 0.5*dR;
    C = (D + sqrt(D^2 - B^2))*0.25;
    A = Rmid + C;

    DZ = Rtop - Rmid;
    y = (B - sqrt(B^2 - 8*C*(DZ - 2*C)))/(4*C);
    cc = sqrt(1 - y^2);

    H = -(Zmax*y)/(2*(cc^3));
    G = Zmax*(1 - 2*(y^2))/(cc^3);

    coeffs = [A; B; C; G; H];
end