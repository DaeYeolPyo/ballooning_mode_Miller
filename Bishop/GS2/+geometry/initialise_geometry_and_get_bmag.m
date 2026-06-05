function [bmag, bpol] = initialise_geometry_and_get_bmag(obj, ...
    rhoc, thetaShift, theta)

    [obj.rpmin, obj.rpmax,...
        obj.surf.rmaj, B_T0, avgrmid] = obj.geom.initialise();

    rp = obj.rpofrho(rhoc);
    drhodrp = obj.drho_drp(rp, obj.surf.dr);

    rgrid = obj.rtgrid(rp, theta);
    
end