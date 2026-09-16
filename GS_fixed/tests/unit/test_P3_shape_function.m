clc
clear
close all

%%
tol = 1e-12;

%% ------------------------------------------------------------
% P3 nodal positions in (xi,eta)
% -------------------------------------------------------------

nodes = [
    0,   0;      % 1
    1,   0;      % 2
    0,   1;      % 3

    1/3, 0;      % 4
    2/3, 0;      % 5

    2/3, 1/3;    % 6
    1/3, 2/3;    % 7

    0,   2/3;    % 8
    0,   1/3;    % 9

    1/3, 1/3     % 10
];


%% ------------------------------------------------------------
% Kronecker-delta test
% -------------------------------------------------------------

M = zeros(10,10);

for a = 1:10

    xi  = nodes(a,1);
    eta = nodes(a,2);

    N = construct_P3_shape_functions(xi,eta);

    M(a,:) = N;

end

err = norm(M-eye(10),inf);

fprintf('Kronecker test error = %.3e\n',err);

assert(err < tol, ...
    'P3 Kronecker-delta test failed.');


%% ------------------------------------------------------------
% Partition-of-unity test
% -------------------------------------------------------------

rng(1);

maxErrN = 0;
maxErrD = 0;

for k = 1:1000

    % Generate random point inside reference triangle
    a = rand;
    b = rand;

    if a+b > 1
        a = 1-a;
        b = 1-b;
    end

    xi  = a;
    eta = b;

    [N,dN] = construct_P3_shape_functions(xi,eta);

    maxErrN = max(maxErrN,abs(sum(N)-1));

    maxErrD = max(maxErrD, ...
        max(abs(sum(dN,2))));

end

fprintf('Partition of unity error = %.3e\n',maxErrN);
fprintf('Derivative sum error      = %.3e\n',maxErrD);

assert(maxErrN < tol);
assert(maxErrD < tol);


%% ------------------------------------------------------------
% Cubic interpolation test
% -------------------------------------------------------------
%
% q(xi,eta) is an arbitrary cubic polynomial.

q = @(x,y) ...
      1 ...
    + 2*x ...
    - 3*y ...
    + 4*x.^2 ...
    + 5*x.*y ...
    - 2*y.^2 ...
    + 3*x.^3 ...
    - 4*x.^2.*y ...
    + 1*x.*y.^2 ...
    + 2*y.^3;

qnode = q(nodes(:,1),nodes(:,2));

maxErr = 0;

for k = 1:1000

    a = rand;
    b = rand;

    if a+b > 1
        a = 1-a;
        b = 1-b;
    end

    [N,~] = construct_P3_shape_functions(a,b);

    qh = N*qnode;
    qe = q(a,b);

    maxErr = max(maxErr,abs(qh-qe));

end

fprintf('Cubic reproduction error  = %.3e\n',maxErr);

assert(maxErr < 1e-11, ...
    'Cubic polynomial reproduction failed.');

fprintf('\nAll P3 shape-function tests passed.\n');