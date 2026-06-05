function p = DshapeMillerParams()
    % Miller paper parameters for psiN ~= 0.77.
    p.A = 3.17;
    p.kappa = 1.66;
    p.delta = 0.416;
    p.s_kappa = 0.70;
    p.s_delta = 1.37;
    p.dR0_dr = -0.354;
    p.q = 3.03;
    p.s_hat = 2.47;
    p.alpha = 1.22;

    % Normalization constants. 
    % With r=1 and B0=1, Bp is returned in B0 units.
    p.r = 0.56;
    p.B0 = 1.33;
end
