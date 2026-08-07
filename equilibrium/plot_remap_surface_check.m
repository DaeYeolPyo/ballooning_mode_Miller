function plot_remap_surface_check(surfaces, maps, k)
%PLOT_REMAP_SURFACE_CHECK Visual check for one remapped surface.

    names = ["PEST", "Hamada", "Boozer"];

    figure;

    for n = 1:numel(names)
        name = names(n);

        [thetaRaw, Rraw, Zraw, lambRaw] = get_raw_surface_data(surfaces(k), name);

        map = maps.(name);
        thetaU = map.theta(:);

        Rmap = map.R(k,:).';
        Zmap = map.Z(k,:).';
        lmap = map.lambda(k,:).';

        thetaRaw = thetaRaw(:) - thetaRaw(1);
        thetaRaw = unwrap(thetaRaw(:));

        subplot(3,3,3*n-2);
        plot(thetaRaw, Rraw, 'k.', thetaU, Rmap, 'r-');
        xlabel('\theta'); ylabel('R'); title(name + " R");

        subplot(3,3,3*n-1);
        plot(thetaRaw, Zraw, 'k.', thetaU, Zmap, 'r-');
        xlabel('\theta'); ylabel('Z'); title(name + " Z");

        subplot(3,3,3*n);
        plot(thetaRaw, lambRaw - lambRaw(1), 'k.', thetaU, lmap, 'r-');
        xlabel('\theta'); ylabel('\lambda'); title(name + " lambda");
    end
end