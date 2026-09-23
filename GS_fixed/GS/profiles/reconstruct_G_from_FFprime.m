function profile = reconstruct_G_from_FFprime( ...
    psiN, FFprime, Fboundary)
%RECONSTRUCT_G_FROM_FFPRIME
%
% Reconstruct
%
%   G(psi_N) = F(psi_N)^2
%
% from
%
%   FFprime(psi_N)
%       = F*dF/dpsi_N
%       = 0.5*dG/dpsi_N,
%
% subject to the boundary condition
%
%   F(psi_N = 1) = Fboundary.
%
% Therefore,
%
%   G(psi_N)
%       = Fboundary^2
%         + 2*integral_1^psi_N FFprime(s) ds.
%
% The sign of F is assumed not to change inside the plasma:
%
%   F(psi_N) = sign(Fboundary)*sqrt(G(psi_N)).
%
%
% INPUT
%
%   psiN
%       Normalized-flux profile grid.
%
%       Requirements:
%           - finite
%           - strictly increasing
%           - includes psi_N = 0 and psi_N = 1
%
%       A nonuniform grid is allowed.
%
%   FFprime
%       F*dF/dpsi_N evaluated on psiN.
%
%       FFprime must have the same number of entries as psiN.
%
%   Fboundary
%       Boundary value
%
%           Fboundary = F(psi_N = 1).
%
%       For the current fixed-boundary solver,
%
%           Fboundary = R0*B0.
%
%
% OUTPUT
%
%   profile
%       Structure containing:
%
%       profile.psiN
%           Column vector containing the normalized-flux grid.
%
%       profile.G
%           Reconstructed G = F^2 profile.
%
%       profile.F
%           Signed F profile.
%
%       profile.FFprime
%           Input F*dF/dpsi_N profile.
%
%       profile.dG_dpsiN
%           dG/dpsi_N = 2*FFprime.
%
%       profile.Fboundary
%           Prescribed boundary F.
%
%       profile.Gboundary
%           Prescribed boundary G = Fboundary^2.
%
%       profile.signF
%           Constant sign used to reconstruct F.
%
%       profile.integrated2FFprime
%           Cumulative integral of 2*FFprime from psi_N = 0.
%
%       profile.minimumG
%           Minimum reconstructed G.
%
%       profile.minimumGIndex
%           Index at which minimum G occurs.
%
%       profile.boundaryError
%           Absolute error in G(end) = Fboundary^2.

    narginchk(3, 3);

    %==============================================================
    % Convert profile arrays to columns
    %==============================================================
    psiN = psiN(:);
    FFprime = FFprime(:);

    %==============================================================
    % Input validation
    %==============================================================
    if numel(psiN) ~= numel(FFprime)
        error('GS:GProfile:ProfileSize', ...
            'psiN and FFprime must have equal lengths.');
    end

    if numel(psiN) < 2
        error('GS:GProfile:TooFewSamples', ...
            'At least two profile samples are required.');
    end

    if any(~isfinite(psiN))
        error('GS:GProfile:InvalidPsiN', ...
            'psiN contains NaN or Inf.');
    end

    if any(diff(psiN) <= 0)
        error('GS:GProfile:InvalidPsiNOrdering', ...
            'psiN must be strictly increasing.');
    end

    if any(~isfinite(FFprime))
        error('GS:GProfile:InvalidFFprime', ...
            'FFprime contains NaN or Inf.');
    end

    validateattributes(Fboundary, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonzero'});

    % The current formulation assumes a profile on [0,1].
    endpointTolerance = ...
        100*eps(max(1, max(abs(psiN))));

    if abs(psiN(1)) > endpointTolerance
        error('GS:GProfile:MissingAxisEndpoint', ...
            'psiN must begin at psi_N = 0.');
    end

    if abs(psiN(end)-1) > endpointTolerance
        error('GS:GProfile:MissingBoundaryEndpoint', ...
            'psiN must end at psi_N = 1.');
    end

    % Remove insignificant endpoint roundoff.
    psiN(1) = 0;
    psiN(end) = 1;

    %==============================================================
    % Reconstruct G = F^2
    %
    % cumulativeIntegral(s)
    %     = integral_0^s 2*FFprime(t) dt
    %
    % Therefore,
    %
    % G(s)
    %     = Fboundary^2
    %       + cumulativeIntegral(s)
    %       - cumulativeIntegral(1).
    %==============================================================
    dG_dpsiN = 2*FFprime;

    integrated2FFprime = ...
        cumtrapz(psiN, dG_dpsiN);

    Gboundary = Fboundary^2;

    G = Gboundary ...
        + integrated2FFprime ...
        - integrated2FFprime(end);

    %==============================================================
    % Physical admissibility
    %==============================================================
    if any(~isfinite(G))
        error('GS:GProfile:InvalidG', ...
            'The reconstructed G profile contains NaN or Inf.');
    end

    [minimumG, minimumGIndex] = min(G);

    if minimumG <= 0
        error('GS:GProfile:NonpositiveG', ...
            ['The reconstructed F^2 profile is nonpositive: ', ...
             'min(G) = %.6e at psi_N = %.6g.'], ...
            minimumG, psiN(minimumGIndex));
    end

    %==============================================================
    % Recover signed F
    %
    % This assumes that F does not cross zero in the plasma.
    %==============================================================
    signF = sign(Fboundary);

    F = signF*sqrt(G);

    if any(~isfinite(F))
        error('GS:GProfile:InvalidF', ...
            'The reconstructed F profile contains NaN or Inf.');
    end

    %==============================================================
    % Boundary consistency
    %==============================================================
    boundaryError = abs(G(end)-Gboundary);

    boundaryScale = max([Gboundary, max(abs(G)), 1]);

    if boundaryError > ...
            100*eps(boundaryScale)
        error('GS:GProfile:BoundaryMismatch', ...
            ['Reconstructed G does not satisfy the boundary ', ...
             'condition. Error = %.6e.'], boundaryError);
    end

    % Set the final entries exactly to the prescribed condition.
    G(end) = Gboundary;
    F(end) = Fboundary;

    %==============================================================
    % Package output
    %==============================================================
    profile = struct();

    profile.psiN = psiN;

    profile.G = G;
    profile.F = F;

    profile.FFprime = FFprime;
    profile.dG_dpsiN = dG_dpsiN;

    profile.Fboundary = Fboundary;
    profile.Gboundary = Gboundary;
    profile.signF = signF;

    profile.integrated2FFprime = ...
        integrated2FFprime;

    profile.minimumG = minimumG;
    profile.minimumGIndex = minimumGIndex;

    profile.boundaryError = ...
        abs(G(end)-Gboundary);
end