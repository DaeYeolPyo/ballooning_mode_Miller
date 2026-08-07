function eq_scaled = scale_eq(eq)
    mu0 = 4*pi*1.e-7;

    % Scaling factors
    BN = eq.bcentr;
    aN = eq.rmaxis; % rcentr?
    psi_scale = 1/(BN*(aN^2));
    F_scale = 1/(aN*BN);
    p_scale = mu0/(BN^2);

    eq_scaled = struct();

    % Copy grid points
    eq_scaled.nw = eq.nw;
    eq_scaled.nh = eq.nh;
    eq_scaled.nbbbs = eq.nbbbs;
    eq_scaled.limitr = eq.limitr;

    % Normalize length scales
    eq_scaled.rdim = eq.rdim/aN;
    eq_scaled.zdim = eq.zdim/aN;
    eq_scaled.rcentr = eq.rcentr/aN;
    eq_scaled.rleft = eq.rleft/aN;
    eq_scaled.zmid = eq.zmid/aN;
    eq_scaled.rmaxis = 1.0;
    eq_scaled.zmaxis = eq.zmaxis/aN;
    eq_scaled.rbbbs = eq.rbbbs/aN;
    eq_scaled.zbbbs = eq.zbbbs/aN;
    eq_scaled.rlim = eq.rlim/aN;
    eq_scaled.zlim = eq.zlim/aN;
    eq_scaled.rgrid = eq.rgrid/aN;
    eq_scaled.zgrid = eq.zgrid/aN;

    % Scale flux
    eq_scaled.simag = eq.simag*psi_scale;
    eq_scaled.sibry = eq.sibry*psi_scale;
    eq_scaled.psi_axis = eq.psi_axis*psi_scale;
    eq_scaled.psi_bdry = eq.psi_bdry*psi_scale;
    eq_scaled.psirz = eq.psirz*psi_scale;
    eq_scaled.psin = eq.psin;

    % Normalize magnetic field
    eq_scaled.bcentr = 1.0;

    % Scale equilibrium profiles
    eq_scaled.pres = eq.pres*p_scale;
    eq_scaled.fpol = eq.fpol*F_scale;
    eq_scaled.pprime = eq.pprime*p_scale/psi_scale;
    eq_scaled.ffprim = eq.ffprim*(F_scale^2)/psi_scale;

    % q profile is not scaled
    eq_scaled.qpsi = eq.qpsi;

    % scale Ip
    eq_scaled.current = eq.current*(mu0*BN/aN);
end