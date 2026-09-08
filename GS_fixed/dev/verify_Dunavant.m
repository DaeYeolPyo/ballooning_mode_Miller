clc
clear
close all

% Verification for Gauss_Dunavant_quadrature(p).
%
% Tests Dunavant rules for p = 1,...,20.
%
% Checks:
%
%   1. Number of quadrature points
%   2. Weight normalization
%   3. Barycentric-coordinate consistency
%   4. Polynomial exactness:
%
%          integral_T xi^i eta^j dA
%
%      for all
%
%          i + j <= p
%
%
% Reference triangle:
%
%           eta
%            ^
%            |
%        (0,1)
%          *
%          |\
%          | \
%          |  \
%          |   \
%          *----* ---> xi
%       (0,0) (1,0)
%
%
% Exact monomial integral:
%
%   int_T xi^i eta^j dxi deta
%
%       = i! j! / (i+j+2)!
%
%
% This code automatically detects whether the quadrature routine
% returns:
%
%   sum(w_q) = 1     : normalized Dunavant weights
%
% or
%
%   sum(w_q) = 1/2   : reference-triangle weights.
%
%
% OUTPUT
%
%   result : table containing verification results for p = 1,...,20
%

tol = 1e-12;

% Expected number of points for Dunavant rules p = 1,...,20
expected_nq = [ ...
     1,  3,  4,  6,  7, ...
    12, 13, 16, 19, 25, ...
    27, 33, 37, 42, 48, ...
    52, 61, 70, 73, 79];

% Storage
P                = (1:20).';
Nq               = zeros(20,1);
ExpectedNq       = expected_nq(:);
WeightSum        = zeros(20,1);
MaxMomentError   = zeros(20,1);
Worst_i          = zeros(20,1);
Worst_j          = zeros(20,1);
MinLambda        = zeros(20,1);
NumOutsidePoints = zeros(20,1);
NumNegativeWeight = zeros(20,1);
Passed           = false(20,1);

fprintf('\n');
fprintf('=============================================================\n');
fprintf('       Dunavant Quadrature Verification: p = 1,...,20\n');
fprintf('=============================================================\n\n');

for p = 1:20

    % -------------------------------------------------------------
    % Obtain quadrature rule
    % -------------------------------------------------------------

    [xi_q, eta_q, w_q] = Gauss_Dunavant_quadrature(p);

    xi_q  = xi_q(:);
    eta_q = eta_q(:);
    w_q   = w_q(:);

    nq = length(w_q);

    Nq(p) = nq;

    % Basic dimensional consistency
    if length(xi_q) ~= nq || length(eta_q) ~= nq
        error( ...
            ['p = %d: xi_q, eta_q, and w_q do not have ', ...
             'the same length.'], p);
    end

    % -------------------------------------------------------------
    % Barycentric coordinates
    %
    % lambda1 = 1 - xi - eta
    % lambda2 = xi
    % lambda3 = eta
    % -------------------------------------------------------------

    lambda1 = 1 - xi_q - eta_q;
    lambda2 = xi_q;
    lambda3 = eta_q;

    all_lambda = [lambda1; lambda2; lambda3];

    MinLambda(p) = min(all_lambda);

    % A point is outside the triangle if any barycentric coordinate < 0.
    % Allow tiny negative values caused only by roundoff.
    outside = ...
        (lambda1 < -tol) | ...
        (lambda2 < -tol) | ...
        (lambda3 < -tol);

    NumOutsidePoints(p) = nnz(outside);

    NumNegativeWeight(p) = nnz(w_q < -tol);

    % -------------------------------------------------------------
    % Detect weight convention
    % -------------------------------------------------------------

    sum_w = sum(w_q);

    WeightSum(p) = sum_w;

    if abs(sum_w - 1.0) < 1e-10

        % Dunavant normalized convention:
        %
        %       sum w_q = 1
        %
        % Therefore quadrature approximates
        %
        %       (1 / Area) integral_T f dA
        %
        % Since Area = 1/2,
        %
        % exact normalized moment =
        %
        %       2 * i! j! / (i+j+2)!

        weight_scale = 2.0;

        convention = 'normalized';

    elseif abs(sum_w - 0.5) < 1e-10

        % Reference-triangle convention:
        %
        %       sum w_q = 1/2
        %
        % and weights directly approximate
        %
        %       integral_T f dA

        weight_scale = 1.0;

        convention = 'triangle';

    else

        warning( ...
            'p = %d: unexpected weight sum = %.16e', ...
            p, sum_w);

        % Continue assuming direct triangle integration
        weight_scale = 1.0;

        convention = 'unknown';

    end

    % -------------------------------------------------------------
    % Polynomial exactness test
    %
    % Check every monomial xi^i eta^j with i+j <= p.
    % -------------------------------------------------------------

    max_err = 0.0;
    worstI = 0;
    worstJ = 0;

    for total_degree = 0:p

        for i = 0:total_degree

            j = total_degree - i;

            % Numerical quadrature
            numerical = sum( ...
                w_q .* ...
                xi_q.^i .* ...
                eta_q.^j );

            % Exact integral:
            %
            % integral_T xi^i eta^j dA
            %
            %   = i! j! / (i+j+2)!
            %
            exact_integral = ...
                factorial(i) * factorial(j) ...
                / factorial(i+j+2);

            % Adapt exact value to weight convention
            exact_quadrature_value = ...
                weight_scale * exact_integral;

            err = abs( ...
                numerical - exact_quadrature_value );

            if err > max_err

                max_err = err;

                worstI = i;
                worstJ = j;

            end

        end

    end

    MaxMomentError(p) = max_err;
    Worst_i(p) = worstI;
    Worst_j(p) = worstJ;

    % -------------------------------------------------------------
    % Overall pass/fail
    % -------------------------------------------------------------

    nq_ok = (nq == expected_nq(p));

    if strcmp(convention,'normalized')

        weight_ok = abs(sum_w - 1.0) < tol;

    elseif strcmp(convention,'triangle')

        weight_ok = abs(sum_w - 0.5) < tol;

    else

        weight_ok = false;

    end

    moment_ok = max_err < tol;

    Passed(p) = nq_ok && weight_ok && moment_ok;

    % -------------------------------------------------------------
    % Print result
    % -------------------------------------------------------------

    if Passed(p)
        status = 'PASS';
    else
        status = 'FAIL';
    end

    fprintf( ...
        ['p = %2d | nq = %2d | sum(w) = %.15f | ', ...
         'max err = %.3e at xi^%d eta^%d | %s\n'], ...
        p, nq, sum_w, max_err, worstI, worstJ, status);

end

fprintf('\n');
fprintf('=============================================================\n');

if all(Passed)

    fprintf('All Dunavant rules passed the polynomial exactness test.\n');

else

    fprintf('Some Dunavant rules FAILED.\n');

    failed_p = find(~Passed);

    fprintf('Failed degrees: ');
    fprintf('%d ', failed_p);
    fprintf('\n');

end

fprintf('=============================================================\n\n');

% -----------------------------------------------------------------
% Summary table
% -----------------------------------------------------------------

result = table( ...
    P, ...
    Nq, ...
    ExpectedNq, ...
    WeightSum, ...
    MaxMomentError, ...
    Worst_i, ...
    Worst_j, ...
    MinLambda, ...
    NumOutsidePoints, ...
    NumNegativeWeight, ...
    Passed);

disp(result)