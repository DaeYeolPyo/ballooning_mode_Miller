function Bp = Bpol_from_EFIT(eq, R, Z)
    % Given R and Z, Evaluate Bp(R, Z)
    % Bp = |grad psi|/R

    psiRZ = eq.psirz';
    [dpdR, dpdZ] = gradient(psiRZ, eq.rgrid, eq.zgrid);
    Bpol = hypot(dpdR, dpdZ)./eq.rgrid;

    Bp_func = griddedInterpolant({eq.zgrid, eq.rgrid}, Bpol, ...
        'spline', 'nearest');

    Bp = Bp_func(Z, R);
end
