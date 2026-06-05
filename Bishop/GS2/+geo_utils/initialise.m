function [psi_0_out, psi_a_out,...
            rmaj, B_T0, avgrmid] = initialise(obj)
    psi_a_out = obj.psi_a; psi_0_out = obj.psi_0;
    rmaj = obj.R_mag;
    B_T0 = abs(obj.B_T);
    avgrmid = obj.aminor;

    obj.initialised = true;
end