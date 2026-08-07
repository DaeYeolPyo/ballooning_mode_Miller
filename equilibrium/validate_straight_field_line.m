function val = validate_straight_field_line(eqfunc, surf, ints, ang)
%VALIDATE_SFL_SURFACE Validate one straight-field-line surface.

    R  = surf.R(:);
    Z  = surf.Z(:);
    dl = surf.dl(:);

    out = eqfunc.eval(R, Z);

    gradPsi = out.gradPsi;
    F       = out.F;
    B2      = out.B2;

    q = ang.q_used;

    val = struct();

    % 1. q consistency
    val.q_int    = ints.q_int;
    val.q_geqdsk = ints.q_geqdsk;
    val.q_relerr = abs(ints.q_int - ints.q_geqdsk) / max(1, abs(ints.q_geqdsk));

    % 2. angle monotonicity
    val.theta_PEST_monotone   = all(diff(ang.theta_PEST)   > 0);
    val.theta_Hamada_monotone = all(diff(ang.theta_Hamada) > 0);
    val.theta_Boozer_monotone = all(diff(ang.theta_Boozer) > 0);

    % 3. angle range
    val.theta_PEST_minmax   = [min(ang.theta_PEST),   max(ang.theta_PEST)];
    val.theta_Hamada_minmax = [min(ang.theta_Hamada), max(ang.theta_Hamada)];
    val.theta_Boozer_minmax = [min(ang.theta_Boozer), max(ang.theta_Boozer)];

    % 4. local Jacobian identities on the contour
    JP = local_jacobian(R, gradPsi, ang.theta_PEST, dl);
    JH = local_jacobian(R, gradPsi, ang.theta_Hamada, dl);
    JB = local_jacobian(R, gradPsi, ang.theta_Boozer, dl);

    JP_expected = q * R.^2 ./ F;

    val.PEST_J_relerr = rms(JP - JP_expected) / max(eps, rms(JP_expected));

    val.Hamada_J_relstd = std(JH) / max(eps, abs(mean(JH)));

    JBB2 = JB .* B2;
    val.Boozer_JB2_relstd = std(JBB2) / max(eps, abs(mean(JBB2)));

    % 5. straight-field-line check for Hamada/Boozer
    val.Hamada_qline_relstd = straight_line_q_std( ...
        ang.theta_PEST, ang.theta_Hamada, ang.lambda_Hamada, q);

    val.Boozer_qline_relstd = straight_line_q_std( ...
        ang.theta_PEST, ang.theta_Boozer, ang.lambda_Boozer, q);
end

function J = local_jacobian(R, gradPsi, theta, dl)
% For coordinates (psi, theta, zeta):
%   J = R / |grad psi| * dl/dtheta

    theta = theta(:);
    dl = dl(:);

    dtheta = circshift(theta, -1) - theta;
    dtheta(end) = 2*pi - theta(end) + theta(1);

    dtheta_dl = dtheta ./ dl;

    J = R ./ gradPsi ./ dtheta_dl;
end

function relstd = straight_line_q_std(thetaP, thetaX, lambdaX, q)
% Verify dzeta_X / dtheta_X = q.
%
% Since in PEST:
%   dphi / dthetaP = q
%
% and for X = Hamada or Boozer:
%   zeta_X = phi + lambda_X,
%
% we check:
%   dzeta_X/dtheta_X =
%   (q + dlambda_X/dthetaP) / (dtheta_X/dthetaP)
%   = q.

    thetaP  = thetaP(:);
    thetaX  = thetaX(:);
    lambdaX = lambdaX(:);

    dthetaX_dthetaP = periodic_derivative(thetaX, thetaP);
    dlambda_dthetaP = periodic_derivative(lambdaX, thetaP);

    qline = (q + dlambda_dthetaP) ./ dthetaX_dthetaP;

    relstd = std(qline) / max(eps, abs(mean(qline)));
end

function dydx = periodic_derivative(y, x)
    y = y(:);
    x = x(:);

    yp = circshift(y, -1);
    ym = circshift(y,  1);

    xp = circshift(x, -1);
    xm = circshift(x,  1);

    xp(end) = x(1) + 2*pi;
    xm(1)   = x(end) - 2*pi;

    dydx = (yp - ym) ./ (xp - xm);
end