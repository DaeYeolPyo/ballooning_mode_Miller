clc
clear
close all

% Geometry parameters
R0    = 4.0;
a     = 1.0;
delta = -0.6;
kappa = 1.5;

% Indentation parameters
Ivals = [0.0, 0.2, 0.4, 0.6, 0.8];

% Localization exponent
m = 4;

% Poloidal angle
theta = linspace(0, 2*pi, 1000);

figure;
hold on;

for I = Ivals

    % Conventional Miller NT shape
    thetaM = theta + asin(delta).*sin(theta);

    % Localized outboard indentation
    Find = ((1 + cos(theta))/2).^m;

    % C-shaped boundary
    R = R0 ...
        + a*cos(thetaM) ...
        - I*a.*Find;

    Z = kappa*a*sin(theta);

    plot(R, Z, 'LineWidth', 2, ...
        'DisplayName', sprintf('I = %.1f', I));

end

grid on
box on

xlabel('R [m]')
ylabel('Z [m]')

legend('Location','best')