function sol = solve_eq32_eig(mats, varargin)
%SOLVE_EQ32_EIG Solve Eq. (32) generalized eigenvalue problem
%
%   sol = SOLVE_EQ32_EIG(mats)
%   sol = SOLVE_EQ32_EIG(mats, 'Name', value, ...)
%
% Solves
%   A f = lambda M f
%
% and returns eigenvalues/eigenvectors sorted by real(lambda).
%
% Name-value options
%   'NumModes'    : number of modes for eigs, default 6
%   'UseEigs'     : true/false, default true
%   'Target'      : eigs sigma target, default 'smallestreal'
%
% Output
%   sol.lambda
%   sol.Vi            interior eigenvectors
%   sol.V             full eigenvectors with Dirichlet boundaries
%   sol.li
%   sol.l

    p = inputParser;
    addParameter(p, 'NumModes', 6, @(x)isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'UseEigs', true, @(x)islogical(x) && isscalar(x));
    addParameter(p, 'Target', 'smallestreal');
    parse(p, varargin{:});
    opt = p.Results;

    A = mats.A;
    M = mats.M;

    if opt.UseEigs
        [Vi, D] = eigs(A, M, opt.NumModes, opt.Target);
        lambda = diag(D);
    else
        [Vi, D] = eig(full(A), full(M));
        lambda = diag(D);
    end

    % sort by real part
    [~, idx] = sort(real(lambda), 'ascend');
    lambda = lambda(idx);
    Vi = Vi(:, idx);

    % rebuild full vectors with Dirichlet BC
    n = mats.n;
    nm = size(Vi, 2);
    V = zeros(n, nm);
    V(2:end-1, :) = Vi;

    sol = struct();
    sol.lambda = lambda;
    sol.Vi = Vi;
    sol.V = V;
    sol.li = mats.li;
    sol.l = mats.l;
end