function Bp = Bpol_from_EFIT(eq, R, Z)
    % Given R and Z, Evaluate Bp(R, Z)
    % Bp = |grad psi|/R

    dR = eq.rgrid(2) - eq.rgrid(1);
    dZ = eq.zgrid(2) - eq.zgrid(1);
    [dpdZ, dpdR] = gradient(eq.psirz', dZ, dR);
    Bpol = hypot(dpdR, dpdZ)./eq.rgrid;

    Bp_func = griddedInterpolant({eq.zgrid, eq.rgrid}, Bpol, ...
        'spline', 'nearest');

    Bp = Bp_func(Z, R);
end