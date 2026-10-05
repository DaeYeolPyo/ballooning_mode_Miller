function [lambda, X] = calculate_eigenvalue(thetaDof, Gmat, Cmat, Fmat, bcType, opts)
    arguments
        thetaDof double
        Gmat double
        Cmat double
        Fmat double
        bcType = 'dirichlet'
        opts.NumEigenvalues (1,1) double = Inf
    end
    if isempty(bcType), bcType='dirichlet'; end
    requested=opts.NumEigenvalues;
    if requested<1 || isnan(requested) || (isfinite(requested) && requested~=floor(requested))
        error('SturmLiouville:EigenvalueCount','NumEigenvalues must be a positive integer or Inf.');
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
    A=P.'*Amat*P; B=P.'*Bmat*P;
    scaleMap=speye(size(A,1));
    if isfinite(requested) && requested<size(A,1)-1
        if norm(A-A.','fro')>1.e-12*max(norm(A,'fro'),realmin) || ...
                norm(B-B.','fro')>1.e-12*max(norm(B,'fro'),realmin)
            error('SturmLiouville:Symmetry','Sparse largest-eigenvalue mode requires symmetric matrices.');
        end
        if any(diag(B)<=0)
            error('SturmLiouville:MassMatrix','The mass matrix must have a positive diagonal.');
        end
        scaleMap=spdiags(1./sqrt(diag(B)),0,size(B,1),size(B,1));
        A=scaleMap*A*scaleMap; B=scaleMap*B*scaleMap;
        settings=struct('tol',1.e-10,'maxit',4000,'disp',0, ...
            'v0',ones(size(A,1),1)/sqrt(size(A,1)));
        if requested==1
            % For SPD B, sigma*B-A is positive definite exactly when sigma
            % exceeds the largest generalized eigenvalue. Sparse Cholesky
            % brackets that eigenvalue without confusing a negative mode
            % near zero with a larger positive mode. Shift-invert then
            % recovers its eigenvector efficiently on long theta intervals.
            upper=bracket_largest_eigenvalue(A,B);
            shift=upper+1.e-6*max(1,abs(upper));
            [V,D,flag]=eigs(A,B,1,shift,settings);
        else
            [V,D,flag]=eigs(A,B,requested,'largestreal',settings);
        end
        if flag~=0
            error('SturmLiouville:EigenvalueNonconvergence','Largest eigenvalue iteration did not converge.');
        end
    else
        [V, D] = eig(full(A),full(B));
    end

    lambda = real(diag(D));
    [lambda, order] = sort(real(lambda), 'descend');
    X = P*scaleMap*V(:, order);
end

function upper=bracket_largest_eigenvalue(A,B)
    [~,flag]=chol(B);
    if flag~=0
        error('SturmLiouville:MassMatrix','The mass matrix must be positive definite.');
    end
    [~,flag]=chol(-A);
    step=1.e-6;
    if flag==0
        upper=0; lower=-step;
        for k=1:80
            [~,flag]=chol(lower*B-A);
            if flag~=0, break; end
            lower=2*lower;
        end
        if flag==0, error('SturmLiouville:Bracket','Could not bracket the largest eigenvalue.'); end
    else
        lower=0; upper=step;
        for k=1:80
            [~,flag]=chol(upper*B-A);
            if flag==0, break; end
            upper=2*upper;
        end
        if flag~=0, error('SturmLiouville:Bracket','Could not bracket the largest eigenvalue.'); end
    end
    for k=1:80
        if upper-lower<=1.e-10*max(1,abs(upper)), break; end
        mid=(lower+upper)/2;
        [~,flag]=chol(mid*B-A);
        if flag==0, upper=mid; else, lower=mid; end
    end
end
