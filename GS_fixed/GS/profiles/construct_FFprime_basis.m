function basis = construct_FFprime_basis(psiN, nBasis, enforceZeroEdge)
%CONSTRUCT_FFPRIME_BASIS Low-order basis for direct FF' reconstruction.
%
%   basis = construct_FFprime_basis(psiN,nBasis,enforceZeroEdge)
%
% The default basis is
%
%   b_j(s) = (1-s) P_{j-1}(2s-1),
%
% where P_j is a Legendre polynomial.  Consequently every basis function
% vanishes at the fixed boundary.  The first basis function is exactly
%
%   b_1(s) = 1-s,
%
% so the common model FF'(s)=A(1-s) is represented by one coefficient.

    narginchk(2, 3);

    if nargin < 3 || isempty(enforceZeroEdge)
        enforceZeroEdge = true;
    end

    psiN = psiN(:);

    validateattributes(nBasis, {'numeric'}, ...
        {'real','finite','scalar','integer','positive','<=',12});

    validateattributes(enforceZeroEdge, {'logical','numeric'}, ...
        {'real','finite','scalar'});

    if any(~isfinite(psiN)) || any(psiN < 0) || any(psiN > 1)
        error('GS:FFprimeBasis:InvalidPsiN', ...
            'psiN must be finite and contained in [0,1].');
    end

    x = 2*psiN-1;
    legendreValues = zeros(numel(psiN), nBasis);
    legendreValues(:,1) = 1;

    if nBasis >= 2
        legendreValues(:,2) = x;
    end

    for degree = 2:nBasis-1
        legendreValues(:,degree+1) = ( ...
            (2*degree-1)*x.*legendreValues(:,degree) ...
            -(degree-1)*legendreValues(:,degree-1)) / degree;
    end

    if logical(enforceZeroEdge)
        basis = (1-psiN).*legendreValues;
    else
        basis = legendreValues;
    end
end
