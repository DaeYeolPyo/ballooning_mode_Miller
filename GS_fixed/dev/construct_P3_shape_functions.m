function [N, dNref] = construct_P3_shape_functions(xi, eta)
% [CONSTRUCT_P3_SHAPE_FUNCTIONS]
% P3 Lagrange basis functions on the reference triangle.
%
% Reference triangle:
%
%       lambda1 = 1 - xi - eta
%       lambda2 = xi
%       lambda3 = eta
%
%       xi >= 0
%       eta >= 0
%       xi + eta <= 1
%
%
% Local node ordering:
%
%                   3
%                   *
%                  / \
%                 8   7
%                /     \
%               9  10   6
%              /    *    \
%             *---*---*---*
%             1   4   5   2
%
%
% OUTPUT
%
%   N       : 1 x 10
%
%   dNref   : 2 x 10
%
%       dNref(1,:) = dN/dxi
%       dNref(2,:) = dN/deta

    L1 = 1 - xi - eta;
    L2 = xi;
    L3 = eta;

    %======================================================
    % Basis functions
    %======================================================
    N = zeros(1, 10);

    N(1)  = B1(L1);
    N(2)  = B1(L2);
    N(3)  = B1(L3);
    N(4)  = B2(L1, L2);
    N(5)  = B2(L2, L1);
    N(6)  = B2(L2, L3);
    N(7)  = B2(L3, L2);
    N(8)  = B2(L3, L1);
    N(9)  = B2(L1, L3);
    N(10) = B3(L1, L2, L3);

    %======================================================
    % Derivatives
    %======================================================
    dN_dxi  = zeros(1, 10);
    dN_deta = zeros(1, 10);

    % xi derivatives
    dN_dxi(1)  = -D1(L1);
    dN_dxi(2)  = D1(L2);
    dN_dxi(3)  = 0;
    dN_dxi(4)  = D3(L1) -D2(L1, L2);
    dN_dxi(5)  = D2(L2, L1) - D3(L2);
    dN_dxi(6)  = D2(L2, L3);
    dN_dxi(7)  = D3(L3);
    dN_dxi(8)  = -D3(L3);
    dN_dxi(9)  = -D2(L1, L3);
    dN_dxi(10) = D4(L1, L2, L3);

    % eta derivatives
    dN_deta(1)  = -D1(L1);
    dN_deta(2)  = 0;
    dN_deta(3)  = D1(L3);
    dN_deta(4)  = -D2(L1, L2);
    dN_deta(5)  = -D3(L2);
    dN_deta(6)  = D3(L2);
    dN_deta(7)  = D2(L3, L2);
    dN_deta(8)  = D2(L3, L1) - D3(L3);
    dN_deta(9)  = D3(L1) - D2(L1, L3);
    dN_deta(10) = D4(L1, L3, L2);

    dNref = [dN_dxi; dN_deta];
end

function result = B1(x)
    result = 0.5*x.*(3*x - 1).*(3*x - 2);
end

function result = B2(x, y)
    result = (9/2)*x.*y.*(3*x - 1);
end

function result = B3(x, y, z)
    result = 27*x.*y.*z;
end

function result = D1(x)
    result = (27/2)*x.^2 - 9*x + 1;
end

function result = D2(x, y)
    result = (9/2)*y.*(6*x - 1);
end

function result = D3(x)
    result = (9/2)*x.*(3*x - 1);
end

function result = D4(x, y, z)
    result = 27*z.*(x - y);
end