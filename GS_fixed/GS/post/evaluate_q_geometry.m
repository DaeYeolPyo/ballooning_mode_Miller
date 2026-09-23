function qGeometry = evaluate_q_geometry( ...
    P3, psi, psiAxis, dpsi, ...
    psiNLevels, nR, axisNode, axisMode)

% [EVALUATE_Q_GEOMETRY]
% Evaluate the geometry-dependent part of the
% safety factor:
%   q(psi_N) = F(psi_N)*Qpsi(psi_N),
%   where
%   Qpsi(psi_N) = 1/(2 pi)*oint dl/(R*|grad(psi)|).
%
% This function depends only on the G-S solution and mesh.
% i.e. does not use F, G (=F^2), or q explicitly.
%
% At the magnetic axis, the contour integral is replaced by
%   Qpsi(0) = 1/(R_axis*sqrt(det(H_psi))).
%
% The boundary geometry factor (psi_N=1)
% is not evaluated directly.
% The final boundary q should be extrapolated after multiplying
% the interior geometry factors by F.
%
% <INPUT>
%
%   P3
%       P3 mesh structure.
%
%   psi
%       Ndof x 1 dimensional poloidal-flux vector.
%
%   psiAxis
%       Magnetic-axis flux used to define psi_N.
%
%   dpsi
%       Signed flux span:
%
%           dpsi = psi_boundary - psi_axis.
%
%   psiNLevels
%       Interior normalized-flux levels. All values must satisfy
%
%           0 < psiNLevels < 1.
%
%   nR
%       Number of Cartesian R samples used by contourc.
%
%   axisNode
%       Global P3 node closest to the magnetic axis.
%
%   axisMode
%       'max', 'min', or 'auto'.
%
% <OUTPUT>
%
%   qGeometry
%       Structure containing:
%
%       psiN
%           Completed profile coordinate [0; interior; 1].
%
%       qPerF
%           Geometry factor q/F. The boundary entry is NaN because
%           boundary q is extrapolated after applying F.
%
%       fluxIntegral
%           oint dl/(R*|grad(psi)|). Boundary entry is NaN.
%
%       axis
%           Recovered magnetic-axis information.
%
%       interiorPsiN
%           Interior contour levels.
%
%       interiorQPerF
%           Interior q/F geometry factors.
%
%       minimumGradPsi
%           Minimum |grad(psi)| on each interior contour.
%
%       contourArea
%           Poloidal area enclosed by each selected contour.
%
%       contourPointCount
%           Number of points on each selected contour.
%
%       contourGridSize
%           [nR,nZ].
%
%       valid
%           True if the calculation completed successfully.
    narginchk(8, 8);

    %==============================================================
    % 0. Input validation
    %==============================================================
    requiredMeshFields = { ...
        'points', ...
        'elements', ...
        'P1points', ...
        'P1elements'};

    for k = 1:numel(requiredMeshFields)
        fieldName = requiredMeshFields{k};

        if ~isfield(P3, fieldName)
            error('GS:QGeometry:MissingMeshField', ...
                'P3.%s is required.', fieldName);
        end
    end

    psi = psi(:);

    nDof = size(P3.points, 1);

    if numel(psi) ~= nDof
        error('GS:QGeometry:PsiSize', ...
            'psi must contain one value per P3 node.');
    end

    if any(~isfinite(psi))
        error('GS:QGeometry:InvalidPsi', ...
            'psi contains NaN or Inf.');
    end

    validateattributes(psiAxis, {'numeric'}, ...
        {'real', 'finite', 'scalar'});

    validateattributes(dpsi, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonzero'});

    psiNLevels = psiNLevels(:);

    if numel(psiNLevels) < 3
        error('GS:QGeometry:TooFewSurfaces', ...
            ['At least three interior flux surfaces are required ', ...
             'for boundary-q extrapolation.']);
    end

    if any(~isfinite(psiNLevels)) || ...
            any(psiNLevels <= 0) || ...
            any(psiNLevels >= 1) || ...
            any(diff(psiNLevels) <= 0)
        error('GS:QGeometry:InvalidPsiNLevels', ...
            ['psiNLevels must be finite, strictly increasing, ', ...
             'and contained in the open interval (0,1).']);
    end

    validateattributes(nR, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', '>=', 2});

    validateattributes(axisNode, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'integer', ...
         '>=', 1, '<=', nDof});

    if isstring(axisMode)
        if ~isscalar(axisMode)
            error('GS:QGeometry:AxisMode', ...
                'axisMode must be a scalar string or character vector.');
        end

        axisMode = char(axisMode);
    end

    if ~ischar(axisMode) || ~isrow(axisMode)
        error('GS:QGeometry:AxisMode', ...
            'axisMode must be ''max'', ''min'', or ''auto''.');
    end

    axisMode = lower(strtrim(axisMode));

    if ~ismember(axisMode, {'max', 'min', 'auto'})
        error('GS:QGeometry:AxisMode', ...
            'axisMode must be ''max'', ''min'', or ''auto''.');
    end

    if size(P3.P1points, 2) ~= 2 || ...
            size(P3.P1elements, 2) ~= 3
        error('GS:QGeometry:InvalidP1Mesh', ...
            'P3.P1points and P3.P1elements have incompatible sizes.');
    end

    if any(P3.P1points(:,1) <= 0)
        error('GS:QGeometry:NonpositiveRadius', ...
            'The mesh contains R <= 0.');
    end

    %==============================================================
    % P1 triangulation used for element lookup
    %==============================================================
    P1 = triangulation(P3.P1elements, P3.P1points);

    %==============================================================
    % Magnetic-axis geometry factor
    %
    %   q_axis/F_axis
    %       = 1/(R_axis*sqrt(det(H_psi))).
    %==============================================================
    [qPerFAxis, axisInfo] = evaluate_axis_q_factor( ...
        P1, P3, psi, axisNode, axisMode);

    axisPoint = axisInfo.point;

    %==============================================================
    % Cartesian contour grid
    %==============================================================
    Rlim = [ ...
        min(P3.P1points(:,1)), ...
        max(P3.P1points(:,1))];

    Zlim = [ ...
        min(P3.P1points(:,2)), ...
        max(P3.P1points(:,2))];

    Rspan = diff(Rlim);
    Zspan = diff(Zlim);

    if ~isfinite(Rspan) || Rspan <= 0 || ...
            ~isfinite(Zspan) || Zspan <= 0
        error('GS:QGeometry:InvalidDomainExtent', ...
            'The mesh must have positive R and Z extents.');
    end

    nZ = max(200, round(nR*Zspan/Rspan));

    Rvec = linspace(Rlim(1), Rlim(2), nR);
    Zvec = linspace(Zlim(1), Zlim(2), nZ);

    [Rgrid, Zgrid] = meshgrid(Rvec, Zvec);

    psiGrid = eval_P3_sol( ...
        P1, P3, psi, [Rgrid(:), Zgrid(:)]);

    psiGrid = reshape(psiGrid, size(Rgrid));

    psiNGrid = (psiGrid-psiAxis)/dpsi;

    %==============================================================
    % Interior flux-surface geometry
    %==============================================================
    nSurface = numel(psiNLevels);

    fluxIntegral     = nan(nSurface, 1);
    qPerFInterior    = nan(nSurface, 1);
    minimumGradPsi   = nan(nSurface, 1);
    contourArea      = nan(nSurface, 1);
    contourPointCount = nan(nSurface, 1);

    for k = 1:nSurface
        psiNk = psiNLevels(k);

        contourMatrix = contourc( ...
            Rvec, Zvec, psiNGrid, [psiNk, psiNk]);

        [X, contourInfo] = select_axis_contour( ...
            contourMatrix, axisPoint, Rlim, Zlim);

        if size(X,1) < 2
            error('GS:QGeometry:InvalidContour', ...
                'Invalid contour on psiN = %.6g.', psiNk);
        end

        dX = diff(X, 1, 1);
        dl = hypot(dX(:,1), dX(:,2));

        Xmid = 0.5*(X(1:end-1,:) + X(2:end,:));

        if any(Xmid(:,1) <= 0)
            error('GS:QGeometry:NonpositiveContourRadius', ...
                'A contour on psiN = %.6g contains R <= 0.', psiNk);
        end

        gradPsi = evaluate_P3_gradient_local( ...
            P1, P3, psi, Xmid);

        gradMagnitude = hypot( ...
            gradPsi(:,1), gradPsi(:,2));

        if any(~isfinite(gradMagnitude)) || ...
                any(gradMagnitude <= 0)
            error('GS:QGeometry:InvalidGradient', ...
                ['Invalid |grad(psi)| on psiN = %.6g. ', ...
                 'Increase contour resolution or refine the mesh.'], ...
                psiNk);
        end

        integralValue = sum( ...
            dl./(Xmid(:,1).*gradMagnitude));

        if ~isfinite(integralValue) || integralValue <= 0
            error('GS:QGeometry:InvalidIntegral', ...
                ['The flux-surface integral is invalid on ', ...
                 'psiN = %.6g.'], psiNk);
        end

        fluxIntegral(k)      = integralValue;
        qPerFInterior(k)     = integralValue/(2*pi);
        minimumGradPsi(k)    = min(gradMagnitude);
        contourArea(k)       = contourInfo.area;
        contourPointCount(k) = size(X,1);
    end

    %==============================================================
    % Package geometry information
    %
    % The boundary entry is deliberately NaN. Boundary q is
    % extrapolated after F has been applied to the interior factors.
    % This preserves the behavior of the current evaluate_q_profile.
    %==============================================================
    qGeometry = struct();

    qGeometry.psiN = [0; psiNLevels; 1];

    qGeometry.qPerF = [ ...
        qPerFAxis; ...
        qPerFInterior; ...
        NaN];

    qGeometry.fluxIntegral = [ ...
        2*pi*qPerFAxis; ...
        fluxIntegral; ...
        NaN];

    qGeometry.qPerFAxis = qPerFAxis;
    qGeometry.interiorPsiN = psiNLevels;
    qGeometry.interiorQPerF = qPerFInterior;

    qGeometry.axis = axisInfo;

    qGeometry.minimumGradPsi = minimumGradPsi;
    qGeometry.contourArea = contourArea;
    qGeometry.contourPointCount = contourPointCount;

    qGeometry.contourGridSize = [nR, nZ];
    qGeometry.RGrid = Rvec;
    qGeometry.ZGrid = Zvec;

    qGeometry.inputPsiAxis = psiAxis;
    qGeometry.dpsi = dpsi;

    qGeometry.boundaryMethod = ...
        'extrapolate-q-after-applying-F';

    qGeometry.valid = true;
    qGeometry.failureReason = '';
end


function [qPerFAxis, info] = evaluate_axis_q_factor( ...
    P1, P3, psi, axisNode, axisMode)
%EVALUATE_AXIS_Q_FACTOR
%
% Recover a smooth quadratic representation around the magnetic axis
% and calculate
%
%   q_axis/F_axis = 1/(R_axis*sqrt(det(H_psi))).

    axisPointNodal = P3.points(axisNode,:);

    % Elements containing the nodal-axis DOF.
    centralElements = find( ...
        any(P3.elements == axisNode, 2));

    if isempty(centralElements)
        error('GS:QGeometry:AxisPatch', ...
            'Failed to locate elements around axisNode.');
    end

    % Expanded recovery patch.
    centralNodes = unique( ...
        P3.elements(centralElements,:));

    patchElements = find( ...
        any(ismember(P3.elements, centralNodes), 2));

    patchNodes = unique( ...
        P3.elements(patchElements,:));

    point = P3.points(patchNodes,:);
    value = psi(patchNodes);

    displacement = point-axisPointNodal;

    patchScale = max(hypot( ...
        displacement(:,1), displacement(:,2)));

    if ~isfinite(patchScale) || patchScale <= 0
        error('GS:QGeometry:AxisPatch', ...
            'Invalid magnetic-axis recovery patch.');
    end

    % Scaled coordinates improve conditioning.
    x = displacement(:,1)/patchScale;
    y = displacement(:,2)/patchScale;

    % Quadratic model:
    %
    % psi = c1 + c2*x + c3*y
    %       + 0.5*c4*x^2 + c5*x*y + 0.5*c6*y^2.
    A = [ ...
        ones(size(x)), ...
        x, ...
        y, ...
        0.5*x.^2, ...
        x.*y, ...
        0.5*y.^2];

    radius = hypot(x,y);
    weight = 1./(1 + (radius/0.5).^4);

    sqrtWeight = sqrt(weight);

    weightedA = A.*sqrtWeight;
    weightedValue = value.*sqrtWeight;

    if rank(weightedA) < 6
        error('GS:QGeometry:AxisHessianRank', ...
            'The axis patch cannot determine a quadratic Hessian.');
    end

    coefficient = weightedA\weightedValue;

    gradientAtNodalPoint = ...
        [coefficient(2); coefficient(3)]/patchScale;

    Hpsi = [ ...
        coefficient(4), coefficient(5); ...
        coefficient(5), coefficient(6)]/patchScale^2;

    Hpsi = 0.5*(Hpsi+Hpsi.');

    eigenvalue = eig(Hpsi);

    switch axisMode
        case 'max'
            if any(eigenvalue >= 0)
                error('GS:QGeometry:AxisHessianSign', ...
                    ['Recovered Hessian is not negative definite ', ...
                     'at a maximum-type axis.']);
            end

            resolvedAxisMode = 'max';

        case 'min'
            if any(eigenvalue <= 0)
                error('GS:QGeometry:AxisHessianSign', ...
                    ['Recovered Hessian is not positive definite ', ...
                     'at a minimum-type axis.']);
            end

            resolvedAxisMode = 'min';

        case 'auto'
            if all(eigenvalue < 0)
                resolvedAxisMode = 'max';
            elseif all(eigenvalue > 0)
                resolvedAxisMode = 'min';
            else
                error('GS:QGeometry:AxisHessianSign', ...
                    ['Recovered Hessian is indefinite, so auto ', ...
                     'axis-mode detection failed.']);
            end

        otherwise
            error('GS:QGeometry:AxisMode', ...
                'Unknown axis mode.');
    end

    determinantH = det(Hpsi);

    if ~isfinite(determinantH) || determinantH <= 0
        error('GS:QGeometry:AxisHessianDeterminant', ...
            'Invalid recovered Hessian determinant.');
    end

    % Stationary point of the recovered quadratic model.
    axisShift = -Hpsi\gradientAtNodalPoint;
    axisPointRecovered = axisPointNodal + axisShift.';

    if norm(axisShift) > 0.5*patchScale
        error('GS:QGeometry:AxisShift', ...
            ['Recovered axis lies too far from the nodal-axis ', ...
             'patch center.']);
    end

    if isnan(pointLocation(P1, axisPointRecovered))
        error('GS:QGeometry:AxisOutsideMesh', ...
            'Recovered magnetic axis is outside the mesh.');
    end

    if axisPointRecovered(1) <= 0
        error('GS:QGeometry:AxisRadius', ...
            'Recovered magnetic-axis radius is not positive.');
    end

    qPerFAxis = 1/( ...
        axisPointRecovered(1)*sqrt(determinantH));

    if ~isfinite(qPerFAxis) || qPerFAxis <= 0
        error('GS:QGeometry:AxisFactor', ...
            'Invalid magnetic-axis q/F geometry factor.');
    end

    % Evaluate the recovered quadratic psi at the refined axis.
    u = axisShift(1)/patchScale;
    v = axisShift(2)/patchScale;

    axisBasis = [ ...
        1, ...
        u, ...
        v, ...
        0.5*u^2, ...
        u*v, ...
        0.5*v^2];

    psiAxisRecovered = axisBasis*coefficient;

    fittedValue = A*coefficient;

    fitScale = max( ...
        norm(value-mean(value)), ...
        eps(max(1,norm(value))));

    fitResidual = ...
        norm(fittedValue-value)/fitScale;

    info = struct();

    info.point = axisPointRecovered;
    info.nodalPoint = axisPointNodal;
    info.shift = axisShift.';
    info.psi = psiAxisRecovered;

    info.Hessian = Hpsi;
    info.HessianEigenvalues = eigenvalue;
    info.HessianDeterminant = determinantH;

    info.qPerF = qPerFAxis;
    info.axisMode = resolvedAxisMode;

    info.fitResidual = fitResidual;
    info.patchScale = patchScale;
    info.patchNodeCount = numel(patchNodes);
    info.fitConditionNumber = cond(weightedA);
end

function [X, info] = select_axis_contour( ...
    contourMatrix, axisPoint, Rlim, Zlim)
%SELECT_AXIS_CONTOUR
%
% Select the outermost closed contour enclosing the magnetic axis.

    curves = {};
    column = 1;

    while column <= size(contourMatrix,2)
        nPoint = round(contourMatrix(2,column));
        lastColumn = column+nPoint;

        if nPoint >= 3 && lastColumn <= size(contourMatrix,2)
            curves{end+1} = ... %#ok<AGROW>
                contourMatrix(:,column+1:lastColumn).';
        end

        column = lastColumn+1;
    end

    domainScale = max([diff(Rlim), diff(Zlim), 1]);
    closureTolerance = 1.e-8*domainScale;

    candidates = {};
    areas = [];
    closureGaps = [];

    for k = 1:numel(curves)
        curve = curves{k};

        closureGap = norm( ...
            curve(end,:)-curve(1,:));

        if closureGap > closureTolerance
            continue
        end

        if closureGap > 10*eps(domainScale)
            curve(end+1,:) = curve(1,:);
        end

        enclosesAxis = inpolygon( ...
            axisPoint(1), axisPoint(2), ...
            curve(:,1), curve(:,2));

        if enclosesAxis
            candidates{end+1} = curve; %#ok<AGROW>
            areas(end+1) = ... %#ok<AGROW>
                polyarea(curve(:,1), curve(:,2));
            closureGaps(end+1) = closureGap; %#ok<AGROW>
        end
    end

    if isempty(candidates)
        error('GS:QGeometry:NoAxisContour', ...
            ['No closed flux surface enclosing the magnetic axis ', ...
             'was found. Increase nR or change psiNLevels.']);
    end

    if numel(candidates) > 1
        warning('GS:QGeometry:MultipleAxisContours', ...
            ['Multiple contours enclose the magnetic axis. ', ...
             'The outermost contour will be used.']);
    end

    [selectedArea, index] = max(areas);

    X = candidates{index};

    info = struct();
    info.area = selectedArea;
    info.closureGap = closureGaps(index);
    info.candidateCount = numel(candidates);
end


function gradPsi = evaluate_P3_gradient_local( ...
    P1, P3, psi, Xq)
%EVALUATE_P3_GRADIENT_LOCAL
%
% Evaluate the elementwise P3 gradient at Cartesian query points.

    elementID = pointLocation(P1, Xq);

    if any(isnan(elementID))
        error('GS:QGeometry:PointOutsideMesh', ...
            ['Some flux-surface integration points are outside ', ...
             'the mesh.']);
    end

    nPoint = size(Xq,1);
    gradPsi = nan(nPoint,2);

    for k = 1:nPoint
        e = elementID(k);

        lambda = cartesianToBarycentric( ...
            P1, e, Xq(k,:));

        xi  = lambda(2);
        eta = lambda(3);

        vertexIDs = P3.P1elements(e,:);
        Ve = P3.P1points(vertexIDs,:);

        [~, dNphysical] = ...
            eval_P3_element(xi, eta, Ve);

        solutionIDs = P3.elements(e,:);

        gradient = ...
            dNphysical*psi(solutionIDs);

        gradPsi(k,:) = gradient.';
    end
end