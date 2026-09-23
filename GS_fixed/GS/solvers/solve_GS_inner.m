function [state, result] = solve_GS_inner( ...
    problem, profiles, psi0, options)
%[SOLVE_GS_INNER]
%
% Thin wrapper around Picard_iteration_prepared.
%
% It converts the converged Picard solution into a state structure
% suitable for:
%
%   - q-geometry evaluation
%   - outer G-profile iteration
%   - warm starting the next inner GS solve
%
%
% INPUT
%
%   problem
%       Output of prepare_GS_problem.
%
%   profiles
%       Structure containing:
%
%           profiles.pprime
%           profiles.FFprime
%
%   psi0
%       Optional initial psi. If omitted or empty, use
%
%           problem.initial.psi.
%
%   options
%       Optional Picard options. If omitted or empty, use
%
%           problem.options.Picard.
%
%
% OUTPUT
%
%   state
%       Converged inner GS state.
%
%   result
%       Picard iteration result.

    narginchk(2, 4);

    if ~isstruct(problem) || ~isscalar(problem)
        error('GS:InnerSolve:ProblemType', ...
            'problem must be a scalar structure.');
    end

    if ~isstruct(profiles) || ~isscalar(profiles)
        error('GS:InnerSolve:ProfileType', ...
            'profiles must be a scalar structure.');
    end

    if ~isfield(profiles,'pprime')
        error('GS:InnerSolve:MissingProfile', ...
            'profiles.pprime is required.');
    end

    if ~isfield(profiles,'FFprime')
        error('GS:InnerSolve:MissingProfile', ...
            'profiles.FFprime is required.');
    end

    if nargin < 3 || isempty(psi0)
        if ~isfield(problem,'initial') || ...
                ~isfield(problem.initial,'psi')
            error('GS:InnerSolve:MissingInitialPsi', ...
                'problem.initial.psi is required.');
        end

        psi0 = problem.initial.psi;
    end

    if nargin < 4 || isempty(options)
        if ~isfield(problem,'options') || ...
                ~isfield(problem.options,'Picard')
            options = struct();
        else
            options = problem.options.Picard;
        end
    end

    pprime = profiles.pprime(:);
    FFprime = profiles.FFprime(:);

    if numel(pprime) ~= numel(FFprime)
        error('GS:InnerSolve:ProfileSize', ...
            'pprime and FFprime must have equal lengths.');
    end

    %==============================================================
    % Prepared Picard solve
    %==============================================================
    [psi, result] = Picard_iteration_prepared( ...
        problem, ...
        pprime, ...
        FFprime, ...
        psi0, ...
        options);

    %==============================================================
    % Normalized flux
    %==============================================================
    psiAxis = result.psiAxis;
    dpsi = result.dpsi;
    psiBoundary = result.boundaryValue;

    psiN = ...
        (psi-psiAxis)/dpsi;

    if any(~isfinite(psiN))
        error('GS:InnerSolve:InvalidPsiN', ...
            'The normalized-flux solution contains NaN or Inf.');
    end

    boundaryPsiNError = norm( ...
        psiN(result.boundaryNodes)-1, inf);

    %==============================================================
    % Locate the nodal magnetic axis
    %==============================================================
    freeNodes = result.freeNodes;

    switch result.axisMode
        case 'max'
            [~, localAxisIndex] = ...
                max(psi(freeNodes));

        case 'min'
            [~, localAxisIndex] = ...
                min(psi(freeNodes));

        otherwise
            error('GS:InnerSolve:AxisMode', ...
                'The resolved axis mode must be max or min.');
    end

    axisNode = ...
        freeNodes(localAxisIndex);

    axisPointNodal = ...
        problem.mesh.points(axisNode,:);

    nodalAxisFlux = ...
        psi(axisNode);

    axisFluxDifference = ...
        abs(nodalAxisFlux-psiAxis);

    %==============================================================
    % Package inner state
    %==============================================================
    state = struct();

    state.psi = psi;
    state.psiN = psiN;

    state.psiAxis = psiAxis;
    state.dpsi = dpsi;
    state.psiBoundary = psiBoundary;

    state.axisMode = result.axisMode;
    state.axisNode = axisNode;

    state.axisPoint = axisPointNodal;
    state.axisPointNodal = axisPointNodal;
    state.axisFluxNodal = nodalAxisFlux;
    state.axisFluxDifference = axisFluxDifference;

    state.freeNodes = result.freeNodes;
    state.boundaryNodes = result.boundaryNodes;

    state.finalRHS = result.finalRHS;
    state.finalResidual = result.finalResidual;

    state.pprime = pprime;
    state.FFprime = FFprime;

    state.innerConverged = result.converged;
    state.innerIterations = result.iterations;

    state.boundaryPsiNError = ...
        boundaryPsiNError;
end