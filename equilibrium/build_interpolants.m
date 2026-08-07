function eqfunc = build_interpolants(eq)
    % Build interpolation objects from read_geqdsk output.
    %
    % Input:
    %   eq = read_geqdsk(filename)
    %
    % Output:
    %   eqfunc.psi(R, Z)
    %   eqfunc.psiR(R, Z)
    %   eqfunc.psiZ(R, Z)
    %   eqfunc.F(psiN)
    %   eqfunc.q(psiN)
    %   eqfunc.p(psiN)
    %   eqfunc.pprime(psiN)
    %   eqfunc.FFprime(psiN)

    R = eq.rgrid(:);
    Z = eq.zgrid(:);

    psiRZ = eq.psirz;

    if ~isequal(size(psiRZ), [numel(R), numel(Z)])
        error('Expected eq.psirz size [%d %d], got [%d %d].', ...
            numel(R), numel(Z), size(psiRZ, 1), size(psiRZ, 2));
    end

    [psiZ_grid, psiR_grid] = gradient(psiRZ, Z, R);

    eqfunc = struct();

    eqfunc.psi = griddedInterpolant({R, Z}, psiRZ, 'spline', 'none');
    eqfunc.psiR = griddedInterpolant({R, Z}, psiR_grid, 'spline', 'none');
    eqfunc.psiZ = griddedInterpolant({R, Z}, psiZ_grid, 'spline', 'none');

    psiN_profile = linspace(0, 1, eq.nw).';

    eqfunc.psiN_from_psi = @(psi) (psi - eq.simag)./(eq.sibry - eq.simag);
    eqfunc.psi_from_psiN = @(psiN) eq.simag + psiN.*(eq.sibry - eq.simag);

    eqfunc.F = griddedInterpolant(psiN_profile, eq.fpol);
    eqfunc.p = griddedInterpolant(psiN_profile, eq.pres);
    eqfunc.FFprime = griddedInterpolant(psiN_profile, eq.ffprim);
    eqfunc.pprime = griddedInterpolant(psiN_profile, eq.pprime);
    eqfunc.q = griddedInterpolant(psiN_profile, eq.qpsi);

    eqfunc.eval = @(Rq, Zq) eval_eq_at_RZ(eqfunc, Rq, Zq);
end

function out = eval_eq_at_RZ(eqfunc, Rq, Zq)
    out = struct();

    out.R = Rq;
    out.Z = Zq;

    out.psi = eqfunc.psi(Rq, Zq);
    out.psiR = eqfunc.psiR(Rq, Zq);
    out.psiZ = eqfunc.psiZ(Rq, Zq);
    out.psiN = eqfunc.psiN_from_psi(out.psi);

    out.gradPsi = hypot(out.psiR, out.psiZ);

    out.F = eqfunc.F(out.psiN);
    out.q = eqfunc.q(out.psiN);
    out.p = eqfunc.p(out.psiN);
    out.pprime = eqfunc.pprime(out.psiN);
    out.ffprime = eqfunc.FFprime(out.psiN);

    out.Bp = out.gradPsi./Rq;
    out.Bphi = out.F./Rq;
    out.B2 = out.Bp.^2 + out.Bphi.^2;
end