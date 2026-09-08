function [N, dNphys, detJ, Rq, Zq] = eval_P3_element(xi, eta, Ve)
% [EVAL_P3_ELEMENT]
% Evaluate P3 basis on an affine physical triangle.
%
% INPUT
%
%   xi, eta
%       Reference coordinates.
%
%   Xe
%       3 x 2 physical coordinates of triangle vertices:
%
%           [R1 Z1
%            R2 Z2
%            R3 Z3]
%
%
% OUTPUT
%
%   N
%       1 x 10 P3 basis.
%
%   dNphys
%       2 x 10 physical derivatives:
%
%           row 1 : dN/dR
%           row 2 : dN/dZ
%
%   detJ
%       Jacobian determinant.
%
%   Rq,Zq
%       Physical coordinates of evaluation point.

    [N, dNref] = construct_P3_shape_functions(xi, eta);

    %===============================================
    % Affine geometry Jacobian
    %===============================================
    R1 = Ve(1, 1);
    Z1 = Ve(1, 2);

    R2 = Ve(2, 1);
    Z2 = Ve(2, 2);

    R3 = Ve(3, 1);
    Z3 = Ve(3, 2);

    J = [R2-R1, R3-R1; Z2-Z1, Z3-Z1];

    detJ = det(J);

    if detJ <= 0
        error('Invalid triangle: det(J) <= 0.');
    end

    %===============================================
    % Physical derivatives
    %   grad_RZ N = J^{-T} grad_xieta N
    %===============================================
    dNphys = J'\dNref;

    %===============================================
    % Physical point coordinates
    %   Since the geometry is affine:
    %   R = L1 R1 + L2 R2 + L3 R3
    %===============================================

    L1 = 1 - xi - eta;
    L2 = xi;
    L3 = eta;

    Rq = L1*R1 + L2*R2 + L3*R3;
    Zq = L1*Z1 + L2*Z2 + L3*Z3;
end