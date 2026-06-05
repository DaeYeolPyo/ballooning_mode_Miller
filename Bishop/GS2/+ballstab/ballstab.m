classdef ballstab < handle
    properties
        make_salpha bool
        initialised bool

        n_shat integer
        n_beta integer

        shat_min double
        shat_max double
        theta0   double
        beta_mul double
        beta_div double
        diff     double

        shat_arr   (:, 1) double
        beta_arr   (:, 1) double
        dbdrho_arr (:, 1) double

        stability (:, :) integer
    end

    methods
        function obj = ballstab()
            obj.beta_div = 1.0; obj.beta_mul = 1.0;
            obj.diff = 0.0;
            obj.make_salpha = false;
            obj.n_beta = 1; obj.n_shat = 1;
            obj.shat_max = 0.0; obj.shat_min = 0.0;
            obj.theta0 = 0.0;
        end

        init_ballstab(obj);
    end
end