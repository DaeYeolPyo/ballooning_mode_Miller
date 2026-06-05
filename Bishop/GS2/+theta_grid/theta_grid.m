classdef theta_grid < handle
    properties
        gb_to_cv bool

        bmin        double
        bmax        double
        eps_trapped double
        shat        double
        drhodpsi    double
        kxfac       double
        qval        double
        cvdriftknob double
        gbdriftknob double
        surfarea    double
        dvdrhon     double
        rhoc        double

        ntheta       integer
        ntgrid       integer
        nperiod      integer
        nbset        integer
        eqopt_switch integer

        theta (:, 1) double
        theta2 (:, 1) double
        delthet (:, 1) double
        delthet2 (:, 1) double
        bset (:, 1) double
        bmag (:, 1) double
        gradpar (:, 1) double
        itor_over_B (:, 1) double
        IoB (:, 1) double
        gbdrift (:, 1) double
        gbdrift0 (:, 1) double
        cvdrift (:, 1) double
        cvdrift0 (:, 1) double
        cdrift (:, 1) double
        cdrift0 (:, 1) double
        gds2 (:, 1) double
        gds21 (:, 1) double
        gds22 (:, 1) double
        gds23 (:, 1) double
        gds24 (:, 1) double
        gds24_noq (:, 1) double
        grho (:, 1) double
        jacob (:, 1) double
        Rplot (:, 1) double
        Zplot (:, 1) double
        aplot (:, 1) double
        Bpol (:, 1) double
        Rprime (:, 1) double
        Zprime (:, 1) double
        aprime (:, 1) double
    end

    methods
        function obj = theta_grid(eqopt)
            obj.eqopt_switch = eqopt;
        end
    end
end