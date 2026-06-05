function DD = drho_drp(obj, rp, dr)
    rp1 = rp*(1-dr); rp2 = rp*(1+dr);

    % Selected irho = 3 case only.
    % It can be considered to use irho = 1, but not implemented
    % Please refer to geometry.f90 #1303-1311
    rho1 = obj.psifun(rp1);
    rho2 = obj.psifun(rp2);

    DD = (rho2 - rho1)/(rp2 - rp1);
end