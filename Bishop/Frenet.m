function [u, R] = Frenet(eq, psiN)
    
end

function [R, Z] = extractFluxSurface(eq, psiN_target)
    % Reconstruct target psi from input psiN
    dpsi = eq.sibry - eq.simag;
    psi_target = eq.simag + psiN_target*dpsi;

    % Extract the corresponding contour
    C = contourc(eq.rgrid, eq.zgrid, eq.psirz', [psi_target psi_target]);

    contours = parse_contourc_output(C);

    if isempty(contours)
        error('[extractFluxSurface] No contour found', ...
            'No contour found for psiN_target = %.2f', psiN_target);
    end
end

function contours = parse_contourc_output(C)
    contours = struct('level', {}, 'R', {}, 'Z', {});
    k = 1;
    ncol = size(C, 2);

    while k < ncol
        level = C(1, k);
        npts = C(2, k);

        cols = (k+1):(k+npts);
        if cols(end) > ncol
            break;
        end

        R = C(1, cols).';
        Z = C(2, cols).';

        s = struct();
        s.level = level;
        s.R = R;
        s.Z = Z;

        contours(end+1) = s;
        k = k + npts + 1;
    end
end