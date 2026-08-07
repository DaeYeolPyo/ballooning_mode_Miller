function sl = validate_field_line_straightness(eqfunc, surf, ang)
%VALIDATE_FIELD_LINE_STRAIGHTNESS Check if field lines are straight in angle space.
%
% Uses the geometric field-line advance:
%
%   dphi = F / (R |grad psi|) dl
%
% and checks:
%
%   alpha_X = zeta_X - q theta_X = const.
%
% for PEST, Hamada, and Boozer coordinates.

    R  = surf.R(:);
    Z  = surf.Z(:);
    dl = surf.dl(:);

    out = eqfunc.eval(R, Z);

    gradPsi = out.gradPsi;
    F       = out.F;

    q = ang.q_used;

    phi = cumulative_toroidal_advance(F, R, gradPsi, dl);

    zeta_P = phi;
    zeta_H = phi + ang.lambda_Hamada(:);
    zeta_B = phi + ang.lambda_Boozer(:);

    alpha_P = unwrap(zeta_P - q * ang.theta_PEST(:));
    alpha_H = unwrap(zeta_H - q * ang.theta_Hamada(:));
    alpha_B = unwrap(zeta_B - q * ang.theta_Boozer(:));

    sl = struct();

    sl.psi  = surf.psi;
    sl.psiN = surf.psiN;
    sl.q = q;

    sl.phi = phi;

    sl.zeta_PEST   = zeta_P;
    sl.zeta_Hamada = zeta_H;
    sl.zeta_Boozer = zeta_B;

    sl.alpha_PEST   = alpha_P;
    sl.alpha_Hamada = alpha_H;
    sl.alpha_Boozer = alpha_B;

    sl.PEST_alpha_relstd = std(alpha_P) / max(1, abs(q));
    sl.Hamada_alpha_relstd = std(alpha_H) / max(1, abs(q));
    sl.Boozer_alpha_relstd = std(alpha_B) / max(1, abs(q));

    sl.PEST_linefit_relerr = line_fit_error(ang.theta_PEST(:), zeta_P, q);
    sl.Hamada_linefit_relerr = line_fit_error(ang.theta_Hamada(:), zeta_H, q);
    sl.Boozer_linefit_relerr = line_fit_error(ang.theta_Boozer(:), zeta_B, q);
end

function phi = cumulative_toroidal_advance(F, R, gradPsi, dl)
% Field-line toroidal advance from outboard starting point.
%
% dphi = F / (R |grad psi|) dl

    w = F ./ (R .* gradPsi);

    wp = circshift(w, -1);
    dphi_seg = 0.5 * (w + wp) .* dl;

    phi = [0; cumsum(dphi_seg(1:end-1))];
    phi = unwrap(phi);
end

function relerr = line_fit_error(theta, zeta, q)
% Check zeta = q theta + alpha0.

    theta = theta(:);
    zeta  = unwrap(zeta(:));

    alpha0 = mean(zeta - q * theta);
    zeta_fit = q * theta + alpha0;

    relerr = rms(zeta - zeta_fit) / max(1, rms(zeta));
end