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