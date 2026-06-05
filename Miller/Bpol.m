function Bp = Bpol(eq, psiN, ntheta, coeffs, dpdr)
    if nargin < 4 || isempty(coeffs)
        coeffs = Cshape_metrics(eq, psiN);
    end

    if nargin < 5 || isempty(dpdr)
        dpdr = abs(eq.sibry - eq.simag);
    end

    Bp = Bpol_from_coeffs(coeffs, dpdr, ntheta);
end