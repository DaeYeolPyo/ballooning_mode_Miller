clc
clear
close all

%% 
psiN_list = linspace(0.05, 0.95, 51);

Npoints_list = [256 512 1024 2048];

conv = check_striaght_field_line_convergence('../Miller/geqdsk_circular', psiN_list, Npoints_list);

%%
eq = read_geqdsk('../Miller/geqdsk_circular');
eqfunc = build_interpolants(eq);

target_psiN = 0.5;

surf = extract_flux_surface(eq, target_psiN, Npoints = 1024);
ints = compute_contour_integrals(eq, eqfunc, surf);
ang  = straight_field_line_angles(eqfunc, surf, ints, UseQ = "integral");

sl = validate_field_line_straightness(eqfunc, surf, ang);

fprintf('PEST straightness   = %.3e\n', sl.PEST_linefit_relerr);
fprintf('Hamada straightness = %.3e\n', sl.Hamada_linefit_relerr);
fprintf('Boozer straightness = %.3e\n', sl.Boozer_linefit_relerr);

figure;

plot(ang.theta_PEST, sl.zeta_PEST, 'k.'); hold on
plot(ang.theta_Hamada, sl.zeta_Hamada, 'r.');
plot(ang.theta_Boozer, sl.zeta_Boozer, 'b.');

theta_plot = linspace(0, 2*pi, 200);
plot(theta_plot, sl.q * theta_plot, 'k--');

xlabel('\theta_X');
ylabel('\zeta_X');
legend('PEST', 'Hamada', 'Boozer', 'slope = q', ...
       'Location', 'northwest');
grid on
title('Straight-field-line check in angle space');

figure;

plot(ang.theta_PEST, sl.alpha_PEST - mean(sl.alpha_PEST), 'k-'); hold on
plot(ang.theta_Hamada, sl.alpha_Hamada - mean(sl.alpha_Hamada), 'r-');
plot(ang.theta_Boozer, sl.alpha_Boozer - mean(sl.alpha_Boozer), 'b-');

xlabel('\theta_X');
ylabel('\alpha_X - <\alpha_X>');
legend('PEST', 'Hamada', 'Boozer');
grid on
title('\alpha_X = \zeta_X - q\theta_X constancy');