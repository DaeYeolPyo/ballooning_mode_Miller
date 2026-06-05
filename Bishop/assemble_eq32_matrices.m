function mats = assemble_eq32_matrices(geom, eq29, eq31, varargin)
%ASSEMBLE_EQ32_MATRICES Assemble generalized eigenvalue matrices for Eq. (32)
%
%   mats = ASSEMBLE_EQ32_MATRICES(geom, eq29, eq31)
%   mats = ASSEMBLE_EQ32_MATRICES(..., 'Name', value, ...)
%
% Solves equations of the form
%
%   d/dl ( a(l) dF/dl ) + lambda * m(l) * F = 0
%
% with Dirichlet boundary conditions F(lmin)=F(lmax)=0.
%
% The generalized eigenvalue problem is assembled as
%
%   A * f = lambda * M * f
%
% where
%   A  ~ - d/dl ( a d/dl )
%   M  ~ m(l)
%
% Inputs
%   geom : output from build_bishop_geometry
%   eq29 : output from compute_eq29
%   eq31 : output from compute_eq31
%
% Name-value options
%   'BC'          : 'dirichlet' only for now
%   'UseSparse'   : true/false, default true
%
% Output fields
%   mats.l
%   mats.li              interior grid
%   mats.a               coefficient a(l)
%   mats.m               coefficient m(l)
%   mats.A               stiffness matrix
%   mats.M               mass-like matrix
%   mats.n
%   mats.ni

    p = inputParser;
    addParameter(p, 'BC', 'dirichlet', @(x)ischar(x) || isstring(x));
    addParameter(p, 'UseSparse', true, @(x)islogical(x) && isscalar(x));
    parse(p, varargin{:});
    opt = p.Results;

    bc = lower(string(opt.BC));
    if bc ~= "dirichlet"
        error('assemble_eq32_matrices:UnsupportedBC', ...
            'Only Dirichlet BC is currently supported.');
    end

    %-----------------------------------
    % 1) Pull arrays
    %-----------------------------------
    l = geom.l(:);
    Bp = geom.Bp(:);
    gradS2 = eq29.absGradS2(:);
    B2 = eq31.B2(:);
    curv = eq31.curvatureTerm(:);

    n = numel(l);
    if n < 3
        error('assemble_eq32_matrices:TooFewPoints', ...
            'Need at least 3 points.');
    end

    %-----------------------------------
    % 2) Coefficients for Eq. (32)
    %-----------------------------------
    % a(l) in divergence term
    a = (Bp.^2) .* gradS2 ./ B2;

    % m(l) multiplying eigenvalue lambda = 2*mu0*p'
    m = curv ./ B2;

    % interior points for Dirichlet BC
    li = l(2:end-1);
    ni = numel(li);

    if opt.UseSparse
        A = spalloc(ni, ni, 3*ni);
        M = spalloc(ni, ni, ni);
    else
        A = zeros(ni, ni);
        M = zeros(ni, ni);
    end

    %-----------------------------------
    % 3) Conservative finite-difference assembly
    %-----------------------------------
    for k = 1:ni
        i = k + 1;  % interior point index in full grid

        dlm = l(i)   - l(i-1);
        dlp = l(i+1) - l(i);

        % midpoint coefficient
        a_imh = 0.5 * (a(i-1) + a(i));
        a_iph = 0.5 * (a(i) + a(i+1));

        % row for -d/dl(a dF/dl)
        A(k,k) =  a_imh/dlm + a_iph/dlp;

        if k > 1
            A(k,k-1) = -a_imh/dlm;
        end
        if k < ni
            A(k,k+1) = -a_iph/dlp;
        end

        % generalized eigenvalue weight
        M(k,k) = m(i);
    end

    mats = struct();
    mats.l = l;
    mats.li = li;
    mats.a = a;
    mats.m = m;
    mats.A = A;
    mats.M = M;
    mats.n = n;
    mats.ni = ni;
end