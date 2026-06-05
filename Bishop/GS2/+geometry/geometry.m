classdef geometry < handle
    properties
        irho   integer
        bishop integer
        nth    integer
        ntgrid integer
        big    integer

        rpmin double
        rpmax double

        efit_eq bool

        geom    geo_utils

        surf    struct
        % Corresponding to eikcoefs_output_type
        outputs struct
    end

    methods
        function obj = geometry()
            obj.big = 8;
            obj.efit_eq = true;
            obj.surf = struct();
            obj.outputs = struct();
        end

        ntheta_out = init_uniform_theta_grid(obj, ntheta, nperiod);
        [bmag, bpol] = initialise_geometry_and_get_bmag(obj, ...
            rhoc, thetaShift, theta);
        rp = rpofrho(obj, rho);
        psi = psifun(obj, rp);
        DD = drho_drp(obj, rp, dr);
        rgrid = rtgrid(obj, rp, theta);
    end

    methods (Static)
        [soln, ier] = root(f, fval, a, b, xerrbi, xerrsec);
    end
end