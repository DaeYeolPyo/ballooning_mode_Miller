classdef geo_utilis < handle
    properties
        type_name string

        initialised          bool
        has_full_theta_range bool

        % Number of radial/theta grid points
        nr integer
        nt integer

        psi_0  double
        psi_a  double
        B_T    double
        beta_0 double
        R_mag  double
        Z_mag  double
        aminor double

        % 2D map of coordinates to R, Z, or B
        R_psi (:, :) double
        Z_psi (:, :) double
        B_psi (:, :) double

        % 2D maps of theta and eqpsi
        eqth     (:, :) double
        eqpsi_2d (:, :) double

        % Minor radius gradient in (R, Z, zeta)
        dpcart (:, :, :) double
        % Theta gradient in (R, Z, zeta)
        dtcart (:, :, :) double
        % B gradient in (R, Z, zeta)
        dbcart (:, :, :) double
        % Minor radius gradient in (rho, l, zeta)
        dpbish (:, :, :) double
        % Theta radius gradient in (rho, l, zeta)
        dtbish (:, :, :) double
        % B gradient in (rho, l, zeta)
        dbbish (:, :, :) double

        eqpsi    (:, 1) double
        pressure (:, 1) double
        psi_bar  (:, 1) double
        fp       (:, 1) double
        qsf      (:, 1) double
        beta     (:, 1) double
        diam     (:, 1) double
        rc       (:, 1) double

        eq eeq
    end

    methods
        function obj = geo_utils(eq)
            obj.initialised = false;
            obj.has_full_theta_range = false;

            obj.eq = eq;
        end

        [psi_0_out, psi_a_out,...
            rmaj, B_T0, avgrmid] = initialise(obj);
        calculate_gradients(obj);
        dfm = derm(obj, f, char);
    end

    methods (Static)
        dfcart = eqdcart(dfm, drm, dzm);
        dbish  = eqdbish(dcart, dpcart);
    end
end