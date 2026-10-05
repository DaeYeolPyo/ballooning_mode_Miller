function psiq = eval_P3_sol(P1, P3, psi, Xq)
% [EVAL_P3_SOL]
% Evaluate P3 FEM solution at Cartesian points.
%
% INPUT
%
%   Xq
%       Nq x 2 array [R,Z]
%
% OUTPUT
%
%   psiq
%       Nq x 1
%
% Points outside the triangulation return NaN.  Keep this compatibility
% wrapper while using the common field/gradient evaluator underneath.

    psiq = eval_P3_field(P1, P3, psi, Xq);
end
