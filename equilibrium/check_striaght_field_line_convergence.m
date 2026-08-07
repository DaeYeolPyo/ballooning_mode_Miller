function conv = check_striaght_field_line_convergence(gfile, psiN_list, Npoints_list)
%CHECK_SFL_CONVERGENCE Check convergence of flux-surface/SFL calculations.

    if nargin < 2 || isempty(psiN_list)
        psiN_list = linspace(0.05, 0.90, 18);
    end

    if nargin < 3 || isempty(Npoints_list)
        Npoints_list = [256 512 1024 2048];
    end

    eq = read_geqdsk(gfile);
    eqfunc  = build_interpolants(eq);

    nPsi = numel(psiN_list);
    nN   = numel(Npoints_list);

    q_int    = nan(nPsi, nN);
    q_geqdsk = nan(nPsi, nN);
    V        = nan(nPsi, nN);
    Vprime   = nan(nPsi, nN);

    pest_err   = nan(nPsi, nN);
    ham_J_std  = nan(nPsi, nN);
    boo_JB_std = nan(nPsi, nN);
    ham_q_std  = nan(nPsi, nN);
    boo_q_std  = nan(nPsi, nN);

    for j = 1:nN
        Np = Npoints_list(j);

        fprintf('\n=== Npoints = %d ===\n', Np);

        for k = 1:nPsi
            psiN = psiN_list(k);

            surf = extract_flux_surface(eq, psiN, Npoints = Np);
            ints = compute_contour_integrals(eq, eqfunc, surf);
            ang  = straight_field_line_angles(eqfunc, surf, ints, UseQ = "integral");
            val  = validate_straight_field_line(eqfunc, surf, ints, ang);

            q_int(k,j)    = ints.q_int;
            q_geqdsk(k,j) = ints.q_geqdsk;
            V(k,j)        = ints.V;
            Vprime(k,j)   = ints.Vprime;

            pest_err(k,j)   = val.PEST_J_relerr;
            ham_J_std(k,j)  = val.Hamada_J_relstd;
            boo_JB_std(k,j) = val.Boozer_JB2_relstd;
            ham_q_std(k,j)  = val.Hamada_qline_relstd;
            boo_q_std(k,j)  = val.Boozer_qline_relstd;
        end
    end

    psi_list = eq.simag + psiN_list(:) * (eq.sibry - eq.simag);

    dV_dpsi_num = nan(size(V));
    Vprime_num_relerr = nan(size(V));

    for j = 1:nN
        dV_dpsi_num(:,j) = gradient(V(:,j), psi_list);

        Vprime_num_relerr(:,j) = abs(Vprime(:,j) - dV_dpsi_num(:,j)) ./ ...
            max(1, abs(Vprime(:,j)));
    end

    % Compare every resolution against the finest resolution.
    jRef = nN;

    q_conv      = relerr_to_ref(q_int,  q_int(:,jRef));
    V_conv      = relerr_to_ref(V,      V(:,jRef));
    Vprime_conv = relerr_to_ref(Vprime, Vprime(:,jRef));

    conv = struct();

    conv.gfile = gfile;
    conv.psiN_list = psiN_list(:);
    conv.psi_list  = psi_list;
    conv.Npoints_list = Npoints_list(:);

    conv.q_int = q_int;
    conv.q_geqdsk = q_geqdsk;
    conv.V = V;
    conv.Vprime = Vprime;
    conv.dV_dpsi_num = dV_dpsi_num;

    conv.metric = struct();
    conv.metric.q_vs_geqdsk = abs(q_int - q_geqdsk) ./ max(1, abs(q_geqdsk));
    conv.metric.Vprime_vs_dVdpsi = Vprime_num_relerr;

    conv.metric.PEST_J_relerr = pest_err;
    conv.metric.Hamada_J_relstd = ham_J_std;
    conv.metric.Boozer_JB2_relstd = boo_JB_std;
    conv.metric.Hamada_qline_relstd = ham_q_std;
    conv.metric.Boozer_qline_relstd = boo_q_std;

    conv.convergence = struct();
    conv.convergence.q_to_finest = q_conv;
    conv.convergence.V_to_finest = V_conv;
    conv.convergence.Vprime_to_finest = Vprime_conv;

    print_convergence_summary(conv);
    plot_convergence_summary(conv);
end

function E = relerr_to_ref(A, ref)
    E = abs(A - ref) ./ max(1, abs(ref));
end

function print_convergence_summary(conv)
    N = conv.Npoints_list(:);

    fprintf('\n\n===== SFL convergence summary =====\n');
    fprintf('%10s %12s %12s %12s %12s %12s %12s\n', ...
        'Npoints', 'q-ref', 'V-ref', 'Vp-ref', ...
        'HamadaJ', 'BoozerJB2', 'Vp-dV');

    for j = 1:numel(N)
        qerr  = median(conv.convergence.q_to_finest(:,j), 'omitnan');
        Verr  = median(conv.convergence.V_to_finest(:,j), 'omitnan');
        Vperr = median(conv.convergence.Vprime_to_finest(:,j), 'omitnan');

        HJ    = median(conv.metric.Hamada_J_relstd(:,j), 'omitnan');
        BJB   = median(conv.metric.Boozer_JB2_relstd(:,j), 'omitnan');
        Vpnum = median(conv.metric.Vprime_vs_dVdpsi(:,j), 'omitnan');

        fprintf('%10d %12.3e %12.3e %12.3e %12.3e %12.3e %12.3e\n', ...
            N(j), qerr, Verr, Vperr, HJ, BJB, Vpnum);
    end
end

function plot_convergence_summary(conv)
    N = conv.Npoints_list(:);

    qerr  = median(conv.convergence.q_to_finest, 1, 'omitnan');
    Verr  = median(conv.convergence.V_to_finest, 1, 'omitnan');
    Vperr = median(conv.convergence.Vprime_to_finest, 1, 'omitnan');

    HJ  = median(conv.metric.Hamada_J_relstd, 1, 'omitnan');
    BJB = median(conv.metric.Boozer_JB2_relstd, 1, 'omitnan');

    figure;

    loglog(N, qerr,  'o-', 'LineWidth', 1.5); hold on
    loglog(N, Verr,  's-', 'LineWidth', 1.5);
    loglog(N, Vperr, '^-', 'LineWidth', 1.5);
    loglog(N, HJ,    'd-', 'LineWidth', 1.5);
    loglog(N, BJB,   'x-', 'LineWidth', 1.5);

    grid on
    xlabel('N contour points');
    ylabel('median relative error / variation');

    legend( ...
        'q vs finest', ...
        'V vs finest', ...
        'V'' vs finest', ...
        'Hamada J relstd', ...
        'Boozer J B^2 relstd', ...
        'Location', 'southwest');

    title('SFL coordinate convergence');
end