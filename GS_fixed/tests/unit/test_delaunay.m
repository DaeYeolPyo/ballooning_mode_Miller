clc
clear
close all

%% mesh parameters
Ntheta = 300;
NR = 65;
NZ = 65;

%% specify plasma boundary
theta = linspace(0, 2*pi, 300);

R0 = 3.0;
a = 1.0;
delta = 0.6;
kappa = 1.5;

Rb = R0 + a*cos(theta + asin(delta)*sin(theta));
Zb = kappa*a*sin(theta);

Rb(end) = [];
Zb(end) = [];

Pboundary = [Rb(:), Zb(:)];

Nb = length(Rb);

C = [(1:Nb)', [2:Nb,1]'];

%% assign internal points

Rgrid = linspace(min(Rb),max(Rb),NR);
Zgrid = linspace(min(Zb),max(Zb),NZ);

[RR,ZZ] = meshgrid(Rgrid,Zgrid);

inside = inpolygon(RR,ZZ,Rb,Zb);

Pint = [RR(inside), ZZ(inside)];

P = [Pboundary; Pint];

%% generate mesh
mesh = generate_P1_mesh(P, C, true);