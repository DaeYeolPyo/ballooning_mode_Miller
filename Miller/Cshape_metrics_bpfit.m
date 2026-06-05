function coeffs = Cshape_metrics_bpfit(eq, psiN, ntheta)
    % EFIT extraction mode: fit the C-shape coefficients from the flux
    % surface, then choose their radial derivatives so the local Jacobian
    % reproduces EFIT Bp as closely as the five-parameter C-shape basis allows.
    coeffs = Cshape_metrics(eq, psiN);

    theta = linspace(0, 2*pi, ntheta);
    [R, Z] = CshapeParam(coeffs(:, 1), ntheta);
    [Rt, Zt] = theta_derv(coeffs, ntheta);

    Bp_efit = Bpol_from_EFIT(eq, R, Z);
    dpdr = abs(eq.sibry - eq.simag);

    [Rr0, Zr0] = r_derv(coeffs, ntheta);
    jac0 = Rr0.*Zt - Rt.*Zr0;
    jacSign = sign(mean(jac0));
    if jacSign == 0
        jacSign = 1;
    end

    jacTarget = jacSign*dpdr*hypot(Rt, Zt)./(R.*Bp_efit);

    basis = [ ...
        Zt(:), ...
        sin(theta(:)).*Zt(:), ...
        cos(2*theta(:)).*Zt(:), ...
        -Rt(:).*cos(theta(:)), ...
        Rt(:).*sin(2*theta(:)) ...
    ];

    % Relative-Jacobian weighting tracks relative Bp error better than an
    % unweighted fit, since Bp is proportional to 1/abs(jacobian).
    weights = 1./max(abs(jacTarget(:)), eps);
    coeffs(:, 2) = (basis.*weights)\(jacTarget(:).*weights);
end
