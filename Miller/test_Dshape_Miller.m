clc
clear
close all

p = DshapeMillerParams();
ntheta = 401;

[Bp, out] = Bpol_Dshape(p, ntheta);
R = out.R;
Z = out.Z;

arc = cumsum([0, hypot(diff(R), diff(Z))]);
arcTheta = 2*pi*arc/arc(end);

fprintf('Miller D-shape parameters\n');
fprintf('A = %.4f, kappa = %.4f, delta = %.4f\n', p.A, p.kappa, p.delta);
fprintf('s_kappa = %.4f, s_delta = %.4f, dR0/dr = %.4f\n', ...
    p.s_kappa, p.s_delta, p.dR0_dr);
fprintf('q target = %.4f, q check = %.4f\n', p.q, out.q_check);
fprintf('dpsi/dr = %.8g in the chosen normalization\n', out.dpdr);

figure;

subplot(1, 2, 1);
plot(R, Z, 'LineWidth', 2);
axis equal;
grid on;
xlabel('R');
ylabel('Z');
title('Miller D-shape');

subplot(1, 2, 2);
plot(arcTheta, Bp*p.B0, 'LineWidth', 2);
grid on;
xlabel('Equal arc-length angle');
ylabel('B_p');
title('Poloidal magnetic field');
