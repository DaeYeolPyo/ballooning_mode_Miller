function [R, Z] = get_flux_surface(eq, psiN_target)
    R0 = eq.rmaxis;
    Z0 = eq.zmaxis;

    Rbnd = eq.rbbbs;
    Zbnd = eq.zbbbs;

    npts = length(Rbnd);

    psi_target = psiN_target*(eq.sibry - eq.simag) + eq.simag;

    Rgrid = eq.rgrid(:);
    Zgrid = eq.zgrid(:);
    psiRZ = eq.psirz;
    if isequal(size(psiRZ), [numel(Zgrid), numel(Rgrid)])
        psiRZ = psiRZ.';
    end

    psi_fun = griddedInterpolant({Rgrid, Zgrid}, psiRZ - psi_target, ...
        'spline', 'none');

    maxIter = 1000;
    tol = 1.e-8;

    R = zeros(1, npts);
    Z = zeros(1, npts);
    for i = 1:npts
        dZdR = (Zbnd(i) - Z0)/(Rbnd(i) - R0);

        Rleft = R0;
        Rright = Rbnd(i);
        Zleft = Z0;
        Zright = Zbnd(i);

        fa = psi_fun(Rleft, Zleft);
        fb = psi_fun(Rright, Zright);

        Rsol = 0;
        Zsol = 0;
        for iter = 1:maxIter
            if abs(dZdR) < 1.0
                Rmid = (Rleft + Rright)*0.5;
                Zmid = dZdR*(Rmid - R0) + Z0;
            else
                Zmid = (Zleft + Zright)*0.5;
                Rmid = (Zmid - Z0)/dZdR + R0;
            end

            fc = psi_fun(Rmid, Zmid);

            if abs(fc) < tol
                Rsol = Rmid;
                Zsol = Zmid;
                break;
            end

            if fa*fc < 0.0
                if abs(dZdR) < 1.0
                    Rright = Rmid;
                else
                    Zright = Zmid;
                end
                fb = fc;
            else
                if abs(dZdR) < 1.0
                    Rleft = Rmid;
                else
                    Zleft = Zmid;
                end
                fa = fc;
            end
        end

        R(i) = Rsol;
        Z(i) = Zsol;
    end
end
