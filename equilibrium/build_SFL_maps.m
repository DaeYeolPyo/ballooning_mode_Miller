function maps = build_SFL_maps(surfaces, Ntheta)
    if nargin < 2
        Ntheta = 256;
    end

    Nrho = numel(surfaces);

    theta = linspace(0, 2*pi, Ntheta + 1);
    theta(end) = [];

    psi    = nan(Nrho, 1);
    psiN   = nan(Nrho, 1);
    q      = nan(Nrho, 1);
    F      = nan(Nrho, 1);
    pprime = nan(Nrho, 1);
    V      = nan(Nrho, 1);
    Vprime = nan(Nrho, 1);

    maps.PEST   = init_map(Nrho, Ntheta, theta);
    maps.Hamada = init_map(Nrho, Ntheta, theta);
    maps.Boozer = init_map(Nrho, Ntheta, theta);

    for k = 1:Nrho
        surf = surfaces(k).surf;
        ints = surfaces(k).ints;
        ang  = surfaces(k).ang;

        psi(k)    = surf.psi;
        psiN(k)   = surf.psiN;
        q(k)      = ang.q_used;
        F(k)      = surf.F;
        pprime(k) = surf.pprime;
        V(k)      = ints.V;
        Vprime(k) = ints.Vprime;
        
        pest = remap_single_surface( ...
            surf.R, surf.Z, zeros(size(surf.R)), ang.theta_PEST, theta);
        ham  = remap_single_surface( ...
            surf.R, surf.Z, ang.lambda_Hamada, ang.theta_Hamada, theta);
        booz = remap_single_surface( ...
            surf.R, surf.Z, ang.lambda_Boozer, ang.theta_Boozer, theta);

        maps.PEST.R(k, :) = pest.R;
        maps.PEST.Z(k, :) = pest.Z;
        maps.PEST.lambda(k, :) = pest.lambda;

        maps.Hamada.R(k,:)      = ham.R;
        maps.Hamada.Z(k,:)      = ham.Z;
        maps.Hamada.lambda(k,:) = ham.lambda;

        maps.Boozer.R(k,:)      = booz.R;
        maps.Boozer.Z(k,:)      = booz.Z;
        maps.Boozer.lambda(k,:) = booz.lambda;
    end

    names = ["PEST", "Hamada", "Boozer"];

    for n = 1:numel(names)
        name = names(n);

        maps.(name).psi = psi;
        maps.(name).psiN = psiN;
        maps.(name).q = q;
        maps.(name).F = F;
        maps.(name).pprime = pprime;
        maps.(name).V = V;
        maps.(name).Vprime = Vprime;
    end
end

function map = init_map(Nrho, Ntheta, theta)
    map = struct();

    map.theta = theta(:).';
    map.R = nan(Nrho, Ntheta);
    map.Z = nan(Nrho, Ntheta);
    map.lambda = nan(Nrho, Ntheta);
end

function out = remap_single_surface(R, Z, lambda, theta_raw, theta_uniform)
    R = R(:);
    Z = Z(:);
    lambda = lambda(:);
    theta_raw = theta_raw(:);

    % Rescale poloidal angle into 0 - 2pi
    theta_raw = theta_raw - theta_raw(1);
    theta_raw = unwrap(theta_raw);

    [theta_raw, idx] = sort(theta_raw);
    R = R(idx);
    Z = Z(idx);
    lambda = lambda(idx);

    % Remove duplicate or non-increasing angle points
    keep = [true; diff(theta_raw) > 1.e-12];
    theta_raw = theta_raw(keep);
    R = R(keep);
    Z = Z(keep);
    lambda = lambda(keep);

    % Force the periodic endpoint
    theta_ext = [theta_raw; 2*pi];
    R_ext = [R; R(1)];
    Z_ext = [Z; Z(1)];

    % lambda is a periodic toroidal shift. Remove any tiny end mismatch.
    lambda = lambda - lambda(1);
    lambda_ext = [lambda; 0];

    out = struct();

    out.R = interp1(theta_ext, R_ext, theta_uniform, 'pchip');
    out.Z = interp1(theta_ext, Z_ext, theta_uniform, 'pchip');
    out.lambda = interp1(theta_ext, lambda_ext, theta_uniform, 'pchip');

    out.R = out.R(:).';
    out.Z = out.Z(:).';
    out.lambda = out.lambda(:).';
end