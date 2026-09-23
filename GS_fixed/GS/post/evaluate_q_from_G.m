function qResult = evaluate_q_from_G( ...
    qGeometry, psiNTable, GTable, signF)
%EVALUATE_Q_FROM_G
%
% Construct q from a geometry factor and G = F^2:
%
%   F(psi_N) = signF*sqrt(G(psi_N)),
%
%   q(psi_N) = F(psi_N)*Qpsi(psi_N),
%
% where Qpsi = q/F is supplied by evaluate_q_geometry.
%
% To preserve the behavior of the existing evaluate_q_profile:
%
%   1. F is first constructed on psiNTable.
%   2. F is linearly interpolated to the requested surfaces.
%   3. Interior q values are calculated.
%   4. Boundary q is extrapolated from interior q.
%
% INPUT
%
%   qGeometry
%       Output of evaluate_q_geometry.
%
%   psiNTable
%       Strictly increasing profile coordinate containing the
%       interval [0,1].
%
%   GTable
%       G = F^2 values on psiNTable.
%
%   signF
%       Sign of F. Must be +1 or -1.
%
% OUTPUT
%
%   qResult
%       Structure containing the completed q, F, and G profiles.

    narginchk(4, 4);

    %==============================================================
    % Input validation
    %==============================================================
    requiredGeometryFields = { ...
        'psiN', ...
        'qPerF', ...
        'qPerFAxis', ...
        'interiorPsiN', ...
        'interiorQPerF', ...
        'axis', ...
        'minimumGradPsi', ...
        'contourGridSize'};

    for k = 1:numel(requiredGeometryFields)
        fieldName = requiredGeometryFields{k};

        if ~isfield(qGeometry, fieldName)
            error('GS:QFromG:MissingGeometryField', ...
                'qGeometry.%s is required.', fieldName);
        end
    end

    if isfield(qGeometry, 'valid') && ...
            ~logical(qGeometry.valid)
        if isfield(qGeometry, 'failureReason')
            failureReason = qGeometry.failureReason;
        else
            failureReason = 'Unknown q-geometry failure.';
        end

        error('GS:QFromG:InvalidGeometry', ...
            'qGeometry is invalid: %s', failureReason);
    end

    psiNTable = psiNTable(:);
    GTable = GTable(:);

    if numel(psiNTable) ~= numel(GTable)
        error('GS:QFromG:ProfileSize', ...
            'psiNTable and GTable must have equal lengths.');
    end

    if numel(psiNTable) < 2
        error('GS:QFromG:ProfileSize', ...
            'At least two profile samples are required.');
    end

    if any(~isfinite(psiNTable)) || ...
            any(diff(psiNTable) <= 0)
        error('GS:QFromG:InvalidPsiNTable', ...
            'psiNTable must be finite and strictly increasing.');
    end

    profileTolerance = ...
        100*eps(max(1,max(abs(psiNTable))));

    if psiNTable(1) > profileTolerance || ...
            psiNTable(end) < 1-profileTolerance
        error('GS:QFromG:ProfileRange', ...
            'psiNTable must cover the normalized interval [0,1].');
    end

    if any(~isfinite(GTable)) || any(GTable <= 0)
        error('GS:QFromG:InvalidG', ...
            'GTable must contain finite positive values.');
    end

    validateattributes(signF, {'numeric'}, ...
        {'real', 'finite', 'scalar'});

    if ~(signF == 1 || signF == -1)
        error('GS:QFromG:InvalidFSign', ...
            'signF must be +1 or -1.');
    end

    psiNInterior = qGeometry.interiorPsiN(:);
    qPerFInterior = qGeometry.interiorQPerF(:);

    if numel(psiNInterior) ~= numel(qPerFInterior)
        error('GS:QFromG:GeometrySize', ...
            ['qGeometry.interiorPsiN and interiorQPerF ', ...
             'must have equal lengths.']);
    end

    if numel(psiNInterior) < 3
        error('GS:QFromG:TooFewSurfaces', ...
            ['At least three interior surfaces are required ', ...
             'for boundary-q extrapolation.']);
    end

    if any(~isfinite(qPerFInterior)) || ...
            any(qPerFInterior <= 0)
        error('GS:QFromG:InvalidGeometryFactor', ...
            'Interior q/F geometry factors must be finite and positive.');
    end

    if ~isfinite(qGeometry.qPerFAxis) || ...
            qGeometry.qPerFAxis <= 0
        error('GS:QFromG:InvalidAxisFactor', ...
            'The magnetic-axis q/F factor must be finite and positive.');
    end

    %==============================================================
    % Construct signed F on the supplied profile grid
    %==============================================================
    FTable = signF*sqrt(GTable);

    % Interpolate F rather than G. This matches the current
    % evaluate_q_profile implementation.
    Faxis = interp1( ...
        psiNTable, FTable, 0, 'linear');

    FInterior = interp1( ...
        psiNTable, FTable, psiNInterior, 'linear');

    Fboundary = interp1( ...
        psiNTable, FTable, 1, 'linear');

    if ~isfinite(Faxis) || ...
            any(~isfinite(FInterior)) || ...
            ~isfinite(Fboundary)
        error('GS:QFromG:FInterpolation', ...
            'Failed to interpolate F onto the q-evaluation surfaces.');
    end

    %==============================================================
    % Axis and interior q
    %==============================================================
    qAxis = ...
        Faxis*qGeometry.qPerFAxis;

    qInterior = ...
        FInterior.*qPerFInterior;

    if ~isfinite(qAxis) || any(~isfinite(qInterior))
        error('GS:QFromG:InvalidQ', ...
            'The calculated axis or interior q contains NaN or Inf.');
    end

    %==============================================================
    % Boundary q
    %
    % Boundary q is extrapolated after F has been applied. This is
    % intentionally consistent with the current solve_GS behavior.
    %==============================================================
    [qBoundary, extrapolationSpread] = ...
        extrapolate_boundary_q( ...
            psiNInterior, qInterior);

    %==============================================================
    % Completed profiles
    %==============================================================
    completedPsiN = [0; psiNInterior; 1];

    completedF = [ ...
        Faxis; ...
        FInterior; ...
        Fboundary];

    completedG = completedF.^2;

    completedQ = [ ...
        qAxis; ...
        qInterior; ...
        qBoundary];

    % The boundary q/F value is an effective extrapolated value,
    % not a directly evaluated boundary geometry integral.
    effectiveBoundaryQPerF = ...
        qBoundary/Fboundary;

    completedQPerF = [ ...
        qGeometry.qPerFAxis; ...
        qPerFInterior; ...
        effectiveBoundaryQPerF];

    %==============================================================
    % Package result
    %==============================================================
    qResult = struct();

    qResult.psiN = completedPsiN;

    qResult.q = completedQ;
    qResult.F = completedF;
    qResult.G = completedG;

    qResult.qPerF = completedQPerF;

    qResult.qAxis = qAxis;
    qResult.qBoundary = qBoundary;

    qResult.FAxis = Faxis;
    qResult.FBoundary = Fboundary;

    qResult.GAxis = Faxis^2;
    qResult.GBoundary = Fboundary^2;

    qResult.interiorPsiN = psiNInterior;
    qResult.interiorQ = qInterior;
    qResult.interiorF = FInterior;
    qResult.interiorG = FInterior.^2;

    qResult.minimumGradPsi = ...
        qGeometry.minimumGradPsi;

    qResult.axis = qGeometry.axis;
    qResult.axis.F = Faxis;
    qResult.axis.G = Faxis^2;
    qResult.axis.q = qAxis;

    qResult.boundaryExtrapolationSpread = ...
        extrapolationSpread;

    qResult.boundaryQPerFIsExtrapolated = true;

    qResult.contourGridSize = ...
        qGeometry.contourGridSize;
end


function [qBoundary, spread] = ...
    extrapolate_boundary_q(psiN, q)
%EXTRAPOLATE_BOUNDARY_Q
%
% Extrapolate q to psi_N = 1 from the outermost interior surfaces.

    psiN = psiN(:);
    q = q(:);

    if numel(psiN) ~= numel(q)
        error('GS:QFromG:BoundaryProfileSize', ...
            'psiN and q must have equal lengths.');
    end

    if numel(q) < 3
        error('GS:QFromG:BoundarySamples', ...
            'At least three interior q values are required.');
    end

    if any(~isfinite(psiN)) || ...
            any(~isfinite(q)) || ...
            any(diff(psiN) <= 0)
        error('GS:QFromG:BoundaryProfile', ...
            ['Boundary-q extrapolation requires finite, ', ...
             'strictly ordered profile samples.']);
    end

    nFit = min(5, numel(q));
    index = numel(q)-nFit+1:numel(q);

    % Shift the coordinate so the target boundary is x = 0.
    x = psiN(index)-1;
    y = q(index);

    linearCoefficient = polyfit(x, y, 1);
    qBoundaryLinear = polyval(linearCoefficient, 0);

    if nFit >= 4
        quadraticCoefficient = polyfit(x, y, 2);
        qBoundary = polyval(quadraticCoefficient, 0);

        spread = abs(qBoundary-qBoundaryLinear);
    else
        qBoundary = qBoundaryLinear;
        spread = NaN;
    end

    if ~isfinite(qBoundary)
        error('GS:QFromG:BoundaryExtrapolation', ...
            'Boundary-q extrapolation failed.');
    end

    if isfinite(spread) && ...
            spread > 0.05*max(abs(qBoundary), eps)
        warning('GS:QBoundary:ExtrapolationSensitive', ...
            ['Linear and quadratic q-boundary extrapolations ', ...
             'differ by %.3e.'], spread);
    end
end