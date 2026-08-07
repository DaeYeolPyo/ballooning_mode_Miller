clear
close all

%%
Ntheta = 1024;
%%
eq = read_geqdsk('../Miller/geqdsk_circular');
eqfunc = build_interpolants(eq);

psiN_list = linspace(0.05, 0.90, 100);

surfaces = repmat(struct(), numel(psiN_list), 1);

for k = 1:numel(psiN_list)
    target_psiN = psiN_list(k);

    surf = extract_flux_surface(eq, target_psiN, Npoints = 1024);
    ints = compute_contour_integrals(eq, eqfunc, surf);
    ang  = straight_field_line_angles(eqfunc, surf, ints, UseQ = "integral");

    surfaces(k).surf = surf;
    surfaces(k).ints = ints;
    surfaces(k).ang  = ang;
end

maps = build_SFL_maps(surfaces, Ntheta);

remapVal = validate_SFL_remap(surfaces, maps);

plot_remap_surface_check(surfaces, maps, 10);