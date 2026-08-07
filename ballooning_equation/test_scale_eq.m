clc
clear
close all

%%
addpath('../equilibrium/');
eq = read_geqdsk('../Miller/geqdsk_circular');
psiN = eq.psin;
eq = scale_eq(eq);

%%
fig = figure();
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile([2 1]);
contour(eq.rgrid, eq.zgrid, eq.psirz.', 'LineWidth', 0.8);
hold on;
plot(eq.rmaxis, eq.zmaxis, 'rx', 'MarkerSize', 9, 'LineWidth', 1.6);
plot(eq.rbbbs, eq.zbbbs, 'k-', 'LineWidth', 1.2);

grid on;
xlabel('R');
ylabel('Z');

nexttile(2);

pprimeN = eq.pprime/max(eq.pprime);
ffprimN = eq.ffprim/max(eq.ffprim);

plot(psiN, pprimeN, 'b-', 'LineWidth', 1.5);
hold on;
plot(psiN, ffprimN, 'r--', 'LineWidth', 1.5);

grid on;
xlabel('\psi_N');
ylabel('normalized value');
legend('p''', 'FF''', 'Location', 'best');

nexttile(3);

plot(psiN, eq.qpsi, 'k-', 'LineWidth', 1.5);

grid on;
xlabel('\psi_N');
ylabel('q');

nexttile(5);

psiN_j = linspace(0.02, 0.98, 80).';

[jphi_HFS, jphi_LFS, RHFS, RLFS] = compute_midplane_jphi(eq, psiN_j, 4*pi*1.e-7);

plot(psiN_j, jphi_HFS, 'b-', 'LineWidth', 1.5);
hold on;
plot(psiN_j, jphi_LFS, 'r--', 'LineWidth', 1.5);

grid on;
xlabel('\psi_N');
ylabel('j_\phi');
legend('HFS', 'LFS', 'Location', 'best');

nexttile(6);

plot(psiN, eq.pres, 'k-', 'LineWidth', 1.5);

grid on;
xlabel('\psi_N');
ylabel('p');

function [jHFS, jLFS, RHFS, RLFS] = compute_midplane_jphi(eq, psiN_list, mu0)
    psiN_grid = linspace(0, 1, eq.nw).';

    pprimeF = griddedInterpolant(psiN_grid, eq.pprime(:), 'pchip', 'nearest');
    ffprimF = griddedInterpolant(psiN_grid, eq.ffprim(:), 'pchip', 'nearest');

    n = numel(psiN_list);

    jHFS = nan(n,1);
    jLFS = nan(n,1);
    RHFS = nan(n,1);
    RLFS = nan(n,1);

    for k = 1:n
        psiN = psiN_list(k);
        psi = eq.simag + psiN * (eq.sibry - eq.simag);

        C = contourc(eq.rgrid(:), eq.zgrid(:), eq.psirz.', [psi psi]);
        seg = choose_axis_contour(C, eq.rmaxis, eq.zmaxis);

        if isempty(seg)
            continue
        end

        R = seg.R(:);

        RHFS(k) = min(R);
        RLFS(k) = max(R);

        pp = pprimeF(psiN);
        ff = ffprimF(psiN);

        % Grad-Shafranov convention:
        %   Delta*psi = -mu0 R^2 pprime - FFprime
        %   j_phi = R pprime + FFprime/(mu0 R)
        jHFS(k) = RHFS(k) * pp + ff / (mu0 * RHFS(k));
        jLFS(k) = RLFS(k) * pp + ff / (mu0 * RLFS(k));
    end
end

function seg = choose_axis_contour(C, Raxis, Zaxis)
    seg = [];

    k = 1;
    bestArea = -inf;

    while k < size(C, 2)
        n = C(2, k);
        cols = k + (1:n);

        R = C(1, cols).';
        Z = C(2, cols).';

        if numel(R) > 1 && hypot(R(1) - R(end), Z(1) - Z(end)) < ...
                1e-10 * max(1, max(hypot(R, Z)))
            R(end) = [];
            Z(end) = [];
        end

        if numel(R) >= 4
            inside = inpolygon(Raxis, Zaxis, R, Z);
            area = polyarea(R, Z);

            if inside && area > bestArea
                bestArea = area;
                seg = struct('R', R, 'Z', Z);
            end
        end

        k = k + n + 1;
    end
end