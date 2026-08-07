clc
clear
close all

addpath('../equilibrium/');
addpath('../sturm_liouville_solver/');

%% Build a normalized PEST equilibrium map
targetPsiN = 0.50;
Ntheta = 256;

eq = read_geqdsk('../Miller/geqdsk_circular');
eq = scale_eq(eq);
eqfunc = build_interpolants(eq);

% The exact target and two or more neighboring surfaces on each side are
% required by calculate_ballooning_coeffs.
psiNList = unique([linspace(0.05, 0.90, 101), targetPsiN]);
surfaces = repmat(struct(), numel(psiNList), 1);

for k = 1:numel(psiNList)
    surf = extract_flux_surface( ...
        eq, eqfunc, psiNList(k), Npoints = int32(1024));

    ints = compute_contour_integrals(eq, eqfunc, surf);
    ang = straight_field_line_angles( ...
        eqfunc, surf, ints, UseQ = "geqdsk");

    surfaces(k).surf = surf;
    surfaces(k).ints = ints;
    surfaces(k).ang = ang;
end

maps = build_SFL_maps(surfaces, Ntheta);
map = maps.PEST;
metrics = compute_metrics(map, eqfunc);

%% Compute g, c, and f only on the requested surface
dtheta = map.theta(2) - map.theta(1);
thetaExt = -4*pi:dtheta:4*pi;

coeff = calculate_ballooning_coeffs( ...
    map, metrics, targetPsiN, ...
    Theta0 = 0, ThetaExt = thetaExt);

fprintf('target psiN       = %.6f\n', coeff.target_psiN);
fprintf('q                 = %.8g\n', coeff.q);
fprintf('dq/dpsi           = %.8g\n', coeff.qprime);
fprintf('pprime            = %.8g\n', coeff.pprime);
fprintf('min(g), min(f)    = %.3e, %.3e\n', ...
    coeff.diagnostics.minG, coeff.diagnostics.minF);
fprintf('B2 metric rel RMS = %.3e\n', ...
    coeff.diagnostics.B2MetricRelativeRms);

figure;
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(coeff.theta, coeff.g, 'LineWidth', 1.5);
grid on
ylabel('g');

nexttile;
plot(coeff.theta, coeff.c, 'LineWidth', 1.5);
grid on
ylabel('c');

nexttile;
plot(coeff.theta, coeff.f, 'LineWidth', 1.5);
grid on
xlabel('\theta');
ylabel('f');

%% Optional handoff to the Sturm-Liouville solver
[Gmat, Cmat, Fmat, thetaDof] = construct_matrix( ...
    coeff.theta, coeff.g, coeff.c, coeff.f);

[lambda, X] = calculate_eigenvalue( ...
    thetaDof, Gmat, Cmat, Fmat, 'dirichlet');

fprintf('largest eigenvalue = %.8g\n', lambda(1));

%% Plot the eigenfunction associated with the largest eigenvalue
mode = X(:, 1);

% Remove the arbitrary complex phase and normalize the peak amplitude.
[~, imax] = max(abs(mode));
mode = mode .* exp(-1i*angle(mode(imax)));
mode = real(mode);
mode = mode ./ max(abs(mode));

figure;
plot(thetaDof, mode, 'b-', 'LineWidth', 1.7);
hold on
yline(0, 'k:');
grid on

xlabel('\theta');
ylabel('X / max|X|');
title(sprintf('Largest-eigenvalue mode, \\lambda = %.6g', lambda(1)));
