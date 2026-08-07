function COCOS = check_COCOS(eq, sigma_RphiZ, e_Bp)
    arguments
        eq struct
        sigma_RphiZ (1, 1) double = 1.0
        e_Bp (1, 1) double = 0.0
    end

    sigma_Ip = sign(eq.current);
    sigma_B0 = sign(eq.bcentr);
    sigma_dpsi = sign(eq.sibry - eq.simag);
    sigma_q = main_sign(eq.qpsi);
    sigma_F = main_sign(eq.fpol);
    sigma_pprime = main_sign(eq.pprime);

    sigma_Bp = sigma_Ip*sigma_dpsi;
    sigma_rhotp = sigma_Ip*sigma_B0*sigma_q;

    check_F = (sigma_F == sigma_B0);
    check_pprime = (sigma_pprime == -sigma_Ip*sigma_Bp);

    if ~(check_F && check_pprime)
        error("COCOS convention sign is not valid.");
    end

    if sigma_Bp == +1 && sigma_rhotp == +1
        base = 1;
    elseif sigma_Bp == -1 && sigma_rhotp == -1
        base = 3;
    elseif sigma_Bp == +1 && sigma_rhotp == -1
        base = 5;
    elseif sigma_Bp == -1 && sigma_rhotp == +1
        base = 7;
    else
        COCOS = NaN;
        return
    end

    if sigma_RphiZ == -1
        base = base + 1;
    end

    if e_Bp == 1
        base = base + 10;
    end

    COCOS = base;
end

function s = main_sign(x)
    x = x(isfinite(x));
    tol = 1.e-10*max(1, max(abs(x)));
    x = x(abs(x) > tol);

    if isempty(x)
        s = NaN;
    else
        s = sign(median(x));
    end
end