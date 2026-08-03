function [lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat, bcType)
    if nargin < 5 || isempty(bcType)
        bcType = 'dirichlet';
    end

    % Generalized eigenvalue problem
    % Ax = lambda Bx
    Amat = Cmat - Gmat;
    Bmat = Fmat;

    ndof = size(Gmat, 1);

    if numel(thetaDof) ~= ndof
        error('thetaDof length mush match the matrix dimension.');
    end

    % Boundary condition
    switch lower(bcType)
        case {'natural', 'neumann'}
            % gX' = 0
            P = speye(ndof);
        case 'dirichlet'
            % X(left) = X(right) = 0
            free = 2:ndof-1;
            P = sparse(free, 1:numel(free), 1, ndof, numel(free));
        case 'natural-dirichlet'
            free = 1:ndof-1;
            P = sparse(free, 1:numel(free), 1, ndof, numel(free));
        case 'periodic'
            % X(left) = X(right)
            P = speye(ndof, ndof - 1);
            P(end, 1) = 1;
        otherwise
            error('Unknown boundary condition: %s', bcType);
    end

    % Eigenvalue solving
    [V, D] = eig(full(P.'*Amat*P), full(P.'*Bmat*P));

    lambda = real(diag(D));
    [lambda, order] = sort(real(lambda), 'descend');
    X = P*V(:, order);
end