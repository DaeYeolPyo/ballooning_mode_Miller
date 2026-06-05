function ntheta_out = init_uniform_theta_grid(obj, ntheta, nperiod)
    obj.nth = floor(ntheta/2);
    ntheta_out = obj.nth*2;
    obj.ntgrid = (2*nperiod - 1)*obj.nth;

    obj.outputs.theta = linspace(-obj.ntgrid, obj.ntgrid, ...
        2*obj.ntgrid + 1).*(pi/obj.nth);
end