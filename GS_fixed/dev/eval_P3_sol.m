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
% Points outside the triangulation return NaN.

    nq = size(Xq,1);

    psiq = nan(nq,1);

    tid = pointLocation(P1,Xq);

    inside = ~isnan(tid);

    for k = find(inside).'
        e = tid(k);

        lambda = cartesianToBarycentric( ...
            P1,e,Xq(k,:));

        xi  = lambda(2);
        eta = lambda(3);

        N = shapeP3(xi,eta);

        ids = P3.elements(e,:);

        psiq(k) = N * psi(ids);
    end
end