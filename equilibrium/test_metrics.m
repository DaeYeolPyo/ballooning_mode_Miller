clc
clear
close all

%%
Ntheta = 256;
%%
eq = read_geqdsk('../Miller/geqdsk_circular');
eqfunc = build_interpolants(eq);

psiN_list = linspace(0.05, 0.90, 100);

surfaces = repmat(struct(), numel(psiN_list), 1);

for k = 1:numel(psiN_list)
    target_psiN = psiN_list(k);

    surf = extract_flux_surface(eq, eqfunc, target_psiN, Npoints = 1024);
    ints = compute_contour_integrals(eq, eqfunc, surf);
    ang  = straight_field_line_angles(eqfunc, surf, ints, UseQ = "integral");

    surfaces(k).surf = surf;
    surfaces(k).ints = ints;
    surfaces(k).ang  = ang;
end

maps = build_SFL_maps(surfaces, Ntheta);

G_P = compute_metrics(maps.PEST, eqfunc);

%% Check if the metric tensor can regenerate B^2 appropriately
k = 50;

B2 = G_P.B2(k, :);
B2_metric = G_P.B2_metric(k, :);

plot(maps.PEST.theta, B2_metric, 'ro--');
hold on
plot(maps.PEST.theta, B2, 'b-', 'LineWidth', 2);

xlabel('\theta', 'FontSize', 14);
ylabel('B^2', 'FontSize', 14);

legend('Smooth B^2', 'B^2 evaluated by metrics');