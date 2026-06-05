function eikcoefs(obj, ntheta_input, nperiod)
    ntheta = obj.init_uniform_theta_grid(ntheta_input, nperiod);
    obj.surf.nt = ntheta;

    thetaShift = 0.0;
    obj.initialise_geometry_and_get_bmag()
end