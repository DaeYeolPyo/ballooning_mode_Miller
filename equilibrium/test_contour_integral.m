clc
clear
close all

%% 
psiN_list = linspace(0.05, 0.95, 51);

eq = read_geqdsk('../Miller/geqdsk_circular');
eqfunc = build_interpolants(eq);

psi_list = eq.simag + psiN_list*(eq.sibry - eq.simag);

V_list      = zeros(size(psiN_list));
Vprime_list = zeros(size(psiN_list));
q_int_list  = zeros(size(psiN_list));
q_g_list    = zeros(size(psiN_list));
q_err_list  = zeros(size(psiN_list));

for k = 1:numel(psiN_list)
    surf = extract_flux_surface(eq, psiN_list(k), Npoints = 1024);
    ints = compute_contour_integrals(eq, eqfunc, surf);

    V_list(k)      = ints.V;
    Vprime_list(k) = ints.Vprime;
    q_int_list(k)  = ints.q_int;
    q_g_list(k)    = ints.q_geqdsk;
    q_err_list(k)  = ints.q_relerr;
end

dV_dpsi_num = gradient(V_list, psi_list);

figure;
plot(psiN_list, Vprime_list, 'k-', ...
     psiN_list, dV_dpsi_num, 'r--');
legend('contour V''', 'numerical dV/d\psi');

figure;
plot(psiN_list, q_int_list, 'k-', ...
     psiN_list, q_g_list, 'r--');
legend('q contour integral', 'q GEQDSK');