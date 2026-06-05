clc
close all
clear

eq = read_geqdsk('../Bishop/geqdsk_PT0.6');

psin0 = 0.8;

bal = eq_ballooning_coefficients(eq, psin0, ...
    'JacMode', 4, ...      % PEST theta
    'NTheta', 256, ...
    'Dpsin', 1e-3);

figure;
plot(bal.theta, bal.g); hold on;
plot(bal.theta, bal.c);
plot(bal.theta, bal.f);
xlabel('\theta');
legend('g','c','f');
grid on;

sol = solve_ballooning_eigenvalue(bal, ...
    'ThetaB', 5*pi, ...
    'N', 501, ...
    'NEigs', 6, ...
    'Plot', true);

sol.lambda
sol.lambda_refined