function [ntheta, theta, ntgrid, nth] = init_uniform_theta_grid(ntheta, theta, nperiod, ntgrid, nth)
    nth = floor(ntheta/2);
    ntheta = nth*2;

    ntgrid = (2*nperiod - 1)*nth;

    theta_len = 2*ntgrid + 1;
    theta = pi*linspace(-ntgrid, ntgrid, theta_len)/nth;
end