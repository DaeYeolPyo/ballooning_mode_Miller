clc
clear
close all

%%
test_delaunay;
P3 = generate_P3_nodes(mesh);

Cp = 4.e5;
CF = 1.e2;
psi_axis = 1.e-2;

%%

% 1. Quadrature
quad = precompute_quadrature(10);

% 2. Global stiffness matrix
K = assemble_K_matrix(P3,quad, true);

% 3. Constant profiles on psi_N grid
nw = 101;

dp_dpsiN  = Cp*ones(nw,1);
FdF_dpsiN = CF*ones(nw,1);

% Fixed value used only for this linear verification
dpsi = -psi_axis;

% For constant profiles, RHS is independent of psiProbe
psiProbe = zeros(size(P3.points,1),1);

% 4. Global source vector
f = assemble_RHS( ...
    P3,psiProbe,quad, ...
    dp_dpsiN,FdF_dpsiN,dpsi);

% 5. Homogeneous fixed-boundary condition
[KII,rhsI,bc] = apply_BC(P3,K,f,0.0);

% 6. Solve and restore full solution
psi = zeros(bc.Ndof,1);

psi(bc.boundaryNodes) = bc.boundaryValues;
psi(bc.freeNodes) = KII\rhsI;

% Boundary error
boundaryError = norm( ...
    psi(bc.boundaryNodes)-bc.boundaryValues,inf);

% Residual of the original, unmodified FEM system
residual = K(bc.freeNodes,:)*psi - f(bc.freeNodes);

relativeResidual = norm(residual) / ...
    max(norm(f(bc.freeNodes)),1);

% Symmetry and positive definiteness
symmetryError = norm(KII-KII.', 'fro')/norm(KII, 'fro');
permutation = symamd(KII);
[~, cholFlag] = chol(KII(permutation, permutation));

fprintf('Boundary error     = %.3e\n',boundaryError);
fprintf('Relative residual  = %.3e\n',relativeResidual);
fprintf('KII symmetry error = %.3e\n',symmetryError);
fprintf('Cholesky flag      = %d\n',cholFlag);

assert(boundaryError < 1e-13);
assert(relativeResidual < 1e-10);
assert(symmetryError < 1e-13);
assert(cholFlag == 0);