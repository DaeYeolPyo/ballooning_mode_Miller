function [outputs, surf] = eikcoefs(ntheta_input, nperiod)
    outputs.theta = init_uniform_theta_grid(ntheta_input, nperiod);
    nth = floor(ntheta/2);
    ntheta = nth*2;
    ntgrid = (2*nperiod - 1)*nth;

    surf.nt = ntheta_input;

    thetaShift = 0.0;

    
end