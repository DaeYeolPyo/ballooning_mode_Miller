clc
clear
close all

A = 1.8;
B = 0.25;
C = 0.16;
G = 0.8;
H = 0.13;

coeffs = [A; B; C; G; H];
ntheta = 301;
[R, Z] = CshapeParam(coeffs, ntheta);

%% Varying B
figure(1);
[R1, Z1] = CshapeParam([A; 0.3; C; G; H], ntheta);
[R2, Z2] = CshapeParam([A; 0.2; C; G; H], ntheta);
[R3, Z3] = CshapeParam([A; 0.5; C; G; H], ntheta);
[R4, Z4] = CshapeParam([A; 0.7; C; G; H], ntheta);

plot(R2, Z2, 'LineWidth', 1.5);
hold on
plot(R1, Z1, 'LineWidth', 1.5);
plot(R3, Z3, 'LineWidth', 1.5);
plot(R4, Z4, 'LineWidth', 1.5);
grid on
legend('B = 0.2', 'B = 0.3', 'B = 0.5', 'B = 0.7');

%% Varying C
figure(2);
[R1, Z1] = CshapeParam([A; B; 0.2; G; H], ntheta);
[R2, Z2] = CshapeParam([A; B; 0.11; G; H], ntheta);

plot(R2, Z2, 'LineWidth', 1.5);
hold on
plot(R, Z, 'LineWidth', 1.5);
plot(R1, Z1, 'LineWidth', 1.5);
grid on
legend('C = 0.11', 'C = 0.16', 'C = 0.2');

%% Varying G
figure(3);
[R1, Z1] = CshapeParam([A; B; C; 0.9; H], ntheta);
[R2, Z2] = CshapeParam([A; B; C; 0.7; H], ntheta);

plot(R2, Z2, 'LineWidth', 1.5);
hold on
plot(R, Z, 'LineWidth', 1.5);
plot(R1, Z1, 'LineWidth', 1.5);
grid on
legend('G = 0.7', 'G = 0.8', 'G = 0.9');

%% Varying B
figure(4);
[R1, Z1] = CshapeParam([A; B; C; G; 0.2], ntheta);
[R2, Z2] = CshapeParam([A; B; C; G; 0.08], ntheta);

plot(R2, Z2, 'LineWidth', 1.5);
hold on
plot(R, Z, 'LineWidth', 1.5);
plot(R1, Z1, 'LineWidth', 1.5);
grid on
legend('H = 0.08', 'H = 0.13', 'H = 0.2');