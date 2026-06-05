function dfm = derm(obj, f, char)
    % dfm(:, :, 1) => derivative w.r.t psi
    % dfm(:, :, 2) => derivative w.r.t theta
    nr = obj.nr; nt = obj.nt;
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

    if obj.has_full_theta_range
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