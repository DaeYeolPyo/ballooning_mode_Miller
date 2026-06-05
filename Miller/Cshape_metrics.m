function coeffs = Cshape_metrics(eq, psiN)
    % coeffs(:, 1, psiN) = [A, B, C, G, H]
    % coeffs(:, 2, psiN) = d/dr[A, B, C, G, H]
    dr = 1.e-4;
    coeffs = zeros(5, 2);

    p = psiN;

    surf = extract_flux_surface(eq, p);
    surf_p = extract_flux_surface(eq, p+dr);
    surf_m = extract_flux_surface(eq, p-dr);
    coefp = fit_Cshape(surf_p);
    coef = fit_Cshape(surf);
    coefm = fit_Cshape(surf_m);

    coeffs(:, 1) = coef;
    coeffs(:, 2) = (coefp - coefm)/(2*dr);
end