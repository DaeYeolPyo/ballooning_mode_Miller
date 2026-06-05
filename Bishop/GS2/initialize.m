function [eq, psi_0_out, psi_a_out, rmaj, B_T0, avgrmid] = initialize(eq)
    eq = calculate_gradients(eq);
    psi_a_out = eq.psi_a;
    psi_0_out = eq.psi_0;
    rmaj = eq.R_mag;
    B_T0 = abs(eq.B_T);
    avgrmid = eq.aminor;
end

function eq = calculate_gradients(eq)
    nr = eq.nr; nt = eq.nt;
    hftr = eq.has_full_theta_range;

    drm = derm(eq.R_psi, nr, nt, hftr, 'E');
    dzm = derm(eq.Z_psi, nr, nt, hftr, 'O');

    dtmp_m = derm(eq.eqpsi_2d, nr, nt, hftr, 'E');
    eq.dpcart = eqdcart(dtmp_m, drm, dzm);
    eq.dpbish = eqdbish(eq.dpcart, eq.dpcart);

    dtmp_m = derm(eq.B_psi, nr, nt, hftr, 'E');
    eq.dbcart = eqdcart(dtmp_m, drm, dzm);
    eq.dbbish = eqdbish(eq.dbcart, eq.dpcart);

    dtmp_m = derm(eq.eqth, nr, nt, hftr, 'T');
    eq.dtcart = eqdcart(dtmp_m, drm, dzm);
    eq.dtbish = eqdbish(eq.dtcart, eq.dpcart);
end

function dfm = derm(f, nr, nt, has_full_theta_range, char)
    % dfm(:, :, 1) => derivative w.r.t psi
    % dfm(:, :, 2) => derivative w.r.t theta
    dfm = zeros(nr, nt, 2);

    i = 1;
    dfm(i, :, 1) = -0.5*(3*f(i, :) - 4*f(i+1, :) + f(i+2, :));

    i = nr;
    dfm(i, :, 1) = 0.5*(3*f(i, :) - 4*f(i-1, :) + f(i-2, :));

    for i = 2:nr-1
        dfm(i, :, 1) = 0.5*(f(i+1, :) - f(i-1, :));
    end

    for j = 2:nt-1
        dfm(:, j, 2) = 0.5*(f(:, j+1) - f(:, j-1));
    end

    if has_full_theta_range
        dfm(:, 1, 2) = 0.5*(f(:, 2) - f(:, nt-1));
        dfm(:, nt, 2) = dfm(:, 1, 2);

        if strcmp(char, 'T')
            dfm(:, 1, 2) = dfm(:, 1, 2) + pi;
            dfm(:, nt, 2) = dfm(:, nt, 2) + pi;
        end
    else
        switch char
            case 'E'
                dfm(:, 1, 2) = 0.0;
                dfm(:, nt, 2) = 0.0;
            case 'O'
                dfm(:, 1, 2) = f(:, 2);
                dfm(:, nt, 2) = -f(:, nt-1);
            case 'T'
                dfm(:, 1, 2) = f(:, 2);
                dfm(:, nt, 2) = pi - f(:, nt-1);
        end
    end
end

function dfcart = eqdcart(dfm, drm, dzm)
    % dfcart(:, :, 1) => derivative w.r.t R
    % dfcart(:, :, 2) => derivative w.r.t Z
    denom = drm(:, :, 1)*dzm(:, :, 2) - drm(:, :, 2)*dzm(:, :, 1);
    dfcart = zeros(size(dfm));

    dfcart(:, :, 1) = dfm(:, :, 1).*dzm(:, :, 2) - dzm(:, :, 1).*dfm(:, :, 2);
    dfcart(:, :, 2) = -dfm(:, :, 1).*drm(:, :, 2) + drm(:, :, 1).*dfm(:, :, 2);

    dfcart(2:end, :, 1) = dfcart(2:end, :, 1)./denom(2:end, :);
    dfcart(2:end, :, 2) = dfcart(2:end, :, 2)./denom(2:end, :);
end

function dbish = eqdbish(dcart, dpcart)
    % dbish(:, :, 1) => (df/dR dpsi/dR + df/dZ dpsi/dZ)/|grad psi|
    % dbish(:, :, 2) => (-df/dR dpsi/dZ + df/dZ dpsi/dR)/|grad psi|
    % denom => |grad psi|
    denom = hypot(dpcart(2:end, :, 1), dpcart(2:end, :, 2));

    dbish = zeros(shape(dcart));
    dbish(:, :, 1) = dcart(:, :, 1).*dpcart(:, :, 1) + dcart(:, :, 2).*dpcart(:, :, 2);
    dbish(:, :, 2) = -dcart(:, :, 1).*dpcart(:, :, 2) + dcart(:, :, 2).*dpcart(:, :, 1);
    dbish(2:end, :, 1) = dbish(2:end, :, 1)./denom;
    dbish(2:end, :, 2) = dbish(2:end, :, 2)./denom;
end