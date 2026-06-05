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