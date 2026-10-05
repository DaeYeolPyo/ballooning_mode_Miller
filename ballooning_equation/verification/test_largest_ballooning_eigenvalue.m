function test_largest_ballooning_eigenvalue
%TEST_LARGEST_BALLOONING_EIGENVALUE Largest means algebraic, not nearest zero.
    setup_indented_cshape_scan;
    theta=linspace(-pi,pi,41).';
    for potential=[-.5 .5 3]
        g=ones(size(theta)); c=potential*g; f=2*g;
        [G,C,F,t]=construct_matrix(theta,g,c,f);
        [allValues,~]=calculate_eigenvalue(t,G,C,F,'dirichlet');
        [largest,X]=calculate_eigenvalue(t,G,C,F,'dirichlet',NumEigenvalues=1);
        analytic=(potential-.25)/2;
        assert(abs(largest-allValues(1))<1.e-9);
        assert(abs(largest-analytic)<1.e-6);
        free=2:numel(t)-1; v=X(free); A=C(free,free)-G(free,free); B=F(free,free);
        assert(norm(A*v-largest*B*v)/(norm(A*v)+norm(largest*B*v))<1.e-8);
    end
    % A positive, nonconstant mass and variable coefficients retain the
    % generalized largest eigenvalue even with several positive modes.
    [G,C,F,t]=construct_matrix(theta,1+.2*cos(theta), ...
        3+.3*sin(theta),2+.1*cos(theta));
    fullValues=calculate_eigenvalue(t,G,C,F,'dirichlet');
    top=calculate_eigenvalue(t,G,C,F,'dirichlet',NumEigenvalues=3);
    assert(norm(top-fullValues(1:3),inf)<1.e-8);
    [G,C,F,t]=construct_matrix(theta,ones(size(theta)),3*ones(size(theta)),2*ones(size(theta)));
    boundaries={'natural','periodic','natural-dirichlet'};
    expected=[1.5 1.5 (3-1/16)/2];
    for k=1:numel(boundaries)
        allValues=calculate_eigenvalue(t,G,C,F,boundaries{k});
        largest=calculate_eigenvalue(t,G,C,F,boundaries{k},NumEigenvalues=1);
        assert(abs(largest-allValues(1))<1.e-9 && abs(largest-expected(k))<1.e-6);
    end
    fprintf('Largest sparse eigenvalues match dense/analytic stable and unstable cases.\n');
end
