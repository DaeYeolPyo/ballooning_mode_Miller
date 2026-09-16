function equilibrium = solve_GS(input)
    %======================================================
    % Set dependencies
    %======================================================
    setup_solve_GS;
    
    %======================================================
    % Check inputs
    %======================================================
    % Grid numbers
    nBoundary        = input.grid.nBoundary;
    NR               = input.grid.NR;
    NZ               = input.grid.NZ;
    quadratureDegree = input.grid.quadratureDegree;
    nProfile         = input.grid.nProfile;

    % Dimensional inputs
    R0     = input.dim.R0; % major radius
    B0     = input.dim.B0; % vacuum toroidal field at R0
    aMinor = input.dim.aMinor;
    
    % Profile inputs
    % To be developed
    % profile_mode = 1: p' and FF'
    %              = 2: p' and <J_par>
    %              = 3: p' and q
    profile_mode = input.profile.mode;
    pprime       = input.profile.pprime;
    FFprime      = input.profile.FFprime;
    pedge        = input.profile.pedge;

    % Plasma boundary inputs
    R_bnd = input.boundary.R;
    Z_bnd = input.boundary.Z;

    psiBoundary = input.boundary.psiBoundary;

    % Picard iteration inputs
    omega             = input.Picard.omega;
    maxIterations     = input.Picard.maxIterations;
    updateTolerance   = input.Picard.updateTolerance;
    residualTolerance = input.Picard.residualTolerance;
    axisMode          = input.Picard.axisMode;
    symmetryTolerance = input.Picard.symmetryTolerance;
    fluxSpanTolerance = input.Picard.fluxSpanTolerance;

    % Outputs
    verbose   = input.output.verbose;
    checkTime = input.output.checkTime;
    showPlot  = input.output.showPlot;

    %% A.1 Generate interior P1 nodes
    Rgrid = linspace(min(R_bnd), max(R_bnd), NR);
    Zgrid = linspace(min(Z_bnd), max(Z_bnd), NZ);
    [RR, ZZ] = meshgrid(Rgrid, Zgrid);

    [insideDomain, onBoundary] = inpolygon(RR, ZZ, R_bnd, Z_bnd);
    strictlyInside = insideDomain & ~onBoundary;

    interiorPoints = [RR(strictlyInside), ZZ(strictlyInside)];

    % Include the geometric center explicitly.
    centerPoint = [R0, 0.0];

    distanceFromCenter = hypot( ...
        interiorPoints(:, 1) - centerPoint(1), ...
        interiorPoints(:, 2) - centerPoint(2));
    coordinateScale = max([max(R_bnd) - min(R_bnd), ...
        max(Z_bnd) - min(Z_bnd), 1]);
    centerTolerance = 1.e-12*coordinateScale;
    interiorPoints(distanceFromCenter < centerTolerance, :) = [];

    plasmaBoundary = [R_bnd(:), Z_bnd(:)];
    P1points = [plasmaBoundary; interiorPoints; centerPoint];
    boundaryConstraints = [(1:nBoundary).', [2:nBoundary, 1].'];

    %% A.2 Generate constrained P1 mesh
    meshTimer = tic;

    meshP1 = generate_P1_mesh(P1points, boundaryConstraints, false);
    triangleCenters = incenter(meshP1);
    triangleInside = inpolygon(triangleCenters(:, 1), ...
                               triangleCenters(:, 2), ...
                               R_bnd, Z_bnd);

    if any(~triangleInside)
        error(['[SOLVE_GS] The constrianed triangulation ' ...
            'contains triangles outside the prescribed LCFS']);
    end

    P1meshTime = toc(meshTimer);

    %% A.3 Upgrade to P3 solution mesh
    P3timer = tic;

    P3 = generate_P3_nodes(meshP1);

    P3meshTime = toc(P3timer);

    Ndof = size(P3.points, 1);
    Nt = size(P3.elements, 1);

    if verbose
        fprintf('\n');
        fprintf('P1 vertices          = %d\n', ...
            size(P3.P1points,1));
        
        fprintf('P1 triangles         = %d\n',Nt);
        
        fprintf('P3 DOFs              = %d\n',Ndof);
        
        fprintf('P3 boundary DOFs     = %d\n', ...
            numel(P3.boundaryNodes));
    end

    if checkTime
        fprintf('P1 mesh time         = %.3f s\n',P1meshTime);
        fprintf('P3 upgrade time      = %.3f s\n',P3meshTime);
    end

    assert(all(P3.points(:,1) > 0), ...
        'The P3 mesh contains R <= 0.');

    %% A.4 Precompute FEM quadrature
    quad = precompute_quadrature(quadratureDegree);

    if verbose
        fprintf('Quadrature degree    = %d\n', ...
            quadratureDegree);
        
        fprintf('Quadrature points    = %d\n', ...
            numel(quad.w));
    end

    %% B.1 Construct initial flux guess
    psiAxisGuess = 0.5;

    distanceToBoundary = inf(Ndof, 1);

    for k = 1:nBoundary
        distanceToPoint = hypot( ...
            P3.points(:, 1) - R_bnd(k), ...
            P3.points(:, 2) - Z_bnd(k));

        distanceToBoundary = min( ...
            distanceToBoundary, ...
            distanceToPoint);
    end

    maximumDistance = max(distanceToBoundary);

    if maximumDistance <= 0
        error('[SOLVE_GS] Failed to construct the initial flux guess.');
    end

    psi0 = psiAxisGuess*distanceToBoundary/maximumDistance;

    % Enforce the homogeneous fixed-boundary value exactly.
    psi0(P3.boundaryNodes) = psiBoundary;

    %% B.2 Load Picard iteration options
    options = struct();

    options.omega             = omega;
    options.maxIterations     = maxIterations;
    options.updateTolerance   = updateTolerance;
    options.residualTolerance = residualTolerance;
    options.axisMode          = axisMode;
    options.symmetryTolerance = symmetryTolerance;
    options.fluxSpanTolerance = fluxSpanTolerance;
    options.verbose           = verbose;

    %% C.1 Solve nonlinear Grad-Shafranov equation
    solveTimer = tic;

    [psi, PicardResult] = Picard_iteration(P3, quad, ...
        pprime, FFprime, ...
        psi0, options);

    solveTime = toc(solveTimer);

    if ~PicardResult.converged
        warning(['The GS solver did not converge.\n' ...
            'Inspect the iteration history' ...
            'before using this equilibrium.']);
    end

    if checkTime
        fprintf('Solver time         = %.3f s\n', solveTime);
    end

    %% C.2 Construct normalized flux
    psiAxis = PicardResult.psiAxis;
    dpsi    = PicardResult.dpsi;
    psiN    = (psi - psiAxis)/dpsi;

    psiNMinimum = min(psiN);
    psiNMaximum = max(psiN);

    %% C.3 Locate the nodal magnetic axis
    freeNodes = PicardResult.freeNodes;

    switch options.axisMode
        case 'max'
            [psiAxisNodal, localAxisIndex] = max(psi(freeNodes));
        case 'min'
            [psiAxisNodal, localAxisIndex] = min(psi(freeNodes));
        otherwise
            error(['[SOLVE_GS] Nodal axis location requires' ...
                'max or min mode.']);
    end

    axisNode  = freeNodes(localAxisIndex);
    axisPoint = P3.points(axisNode, :);
    %% D.1 Fpol and plasma current evaluation
    mu0 = 4*pi*1.e-7;

    psiNProfile = linspace(0, 1, nProfile);
    integratedP = cumtrapz(psiNProfile, pprime);
    pressure = pedge + integratedP - integratedP(end);

    Fvac = R0*B0;
    integratedFF = cumtrapz(psiNProfile, 2*FFprime);
    F2 = Fvac^2 + integratedFF - integratedFF(end);

    if ~isfinite(Fvac) || Fvac == 0
        error('[SOLVE_GS] R0*B0 must be finite and nonzero.');
    end

    if any(~isfinite(F2)) || any(F2 <= 0)
        error('[SOLVE_GS] Reconstructed F^2 must remain positive.');
    end

    Fpol     = sign(Fvac)*sqrt(F2);

    Ip = sum(PicardResult.finalRHS)/mu0;

    % Interior surfaces used for numerical contour integration.
    % Slightly denser sampling is used near the boundary.
    qInteriorPsiN = unique([ ...
        0.05:0.05:0.80, ...
        0.825:0.025:0.95]).';

    qContourNR = 500;

    qResult = evaluate_q_profile( ...
        P3, psi, psiAxis, dpsi, ...
        qInteriorPsiN, qContourNR, ...
        axisNode, options.axisMode, ...
        psiNProfile, Fpol);

    % Interpolate the completed [0,1] q profile onto the common profile grid.
    qProfile = interp1( ...
        qResult.psiN, qResult.q, ...
        psiNProfile, 'pchip');

    Baxis = Fpol(1)/qResult.axis.point(1);

    %% D.2 Evaluate global parameters
    globalResult = evaluate_global_parameters( ...
        P3, quad, psi, ...
        psiAxis, dpsi, ...
        psiNProfile, pressure, ...
        R0, B0, aMinor, Ip);

    %% E.1 Equilibrium result wrap-up
    equilibrium = struct();

    equilibrium.mesh = P3;

    equilibrium.psi = psi;
    equilibrium.psiN = psiN;

    equilibrium.axisNode = axisNode;

    equilibrium.axisPoint = axisPoint;

    % profiles
    equilibrium.profile = struct();

    equilibrium.profile.psiN      = psiNProfile;
    equilibrium.profile.pressure  = pressure;
    equilibrium.profile.dp_dpsiN  = pprime;
    equilibrium.profile.Fpol      = Fpol;
    equilibrium.profile.FdF_dpsiN = FFprime;
    equilibrium.profile.q         = qProfile;

    % global parameters
    equilibrium.global = struct();

    equilibrium.global.Raxis = axisPoint(1);
    equilibrium.global.Zaxis = axisPoint(2);
    equilibrium.global.simag = psiAxis;
    equilibrium.global.sibry = psiBoundary; % Should be zero
    equilibrium.global.Baxis = Baxis;
    equilibrium.global.Ip    = Ip;

    equilibrium.global.volume  = globalResult.volume;   
    equilibrium.global.betaTor = globalResult.betaTor;
    equilibrium.global.betaPol = globalResult.betaPol;
    equilibrium.global.li3     = globalResult.li3;
    equilibrium.global.betaN   = globalResult.betaN;
    
    equilibrium.solverResult = PicardResult;

    %% E.2 GS solve report
    if verbose
        report_GS_result(equilibrium);
    end

    if showPlot
        plot_GS_equilibrium(equilibrium);
    end
end

function qResult = evaluate_q_profile( ...
    P3, psi, psiAxis, dpsi, psiNLevels, ...
    nR, axisNode, axisMode, psiNTable, FTable)

    % Recovered magnetic axis and semi-analytic q_axis
    Faxis = interp1(psiNTable, FTable, 0, 'linear');

    [qAxis, axisInfo] = evaluate_axis_q_hessian( ...
        P3, psi, axisNode, axisMode, Faxis);

    axisPoint = axisInfo.point;

    % Element locations
    P1 = triangulation(P3.P1elements, P3.P1points);

    % Sample Cartesian grid
    Rlim = [min(P3.P1points(:,1)), max(P3.P1points(:,1))];
    Zlim = [min(P3.P1points(:,2)), max(P3.P1points(:,2))];

    Rspan = diff(Rlim);
    Zspan = diff(Zlim);

    nZ = max(200, round(nR*Zspan/Rspan));

    Rvec = linspace(Rlim(1), Rlim(2), nR);
    Zvec = linspace(Zlim(1), Zlim(2), nZ);

    [Rgrid,Zgrid] = meshgrid(Rvec,Zvec);

    psiGrid = eval_P3_sol( ...
        P1, P3, psi, [Rgrid(:),Zgrid(:)]);

    psiGrid = reshape(psiGrid,size(Rgrid));
    psiNGrid = (psiGrid-psiAxis)/dpsi;

    % Interior flux-surface integrals
    nSurface = numel(psiNLevels);

    qInterior = nan(nSurface,1);
    minimumGradPsi = nan(nSurface,1);

    for k = 1:nSurface
        psiNk = psiNLevels(k);

        C = contourc( ...
            Rvec, Zvec, psiNGrid, [psiNk,psiNk]);

        X = select_axis_contour( ...
            C, axisPoint, Rlim, Zlim);

        dX = diff(X,1,1);
        dl = hypot(dX(:,1),dX(:,2));

        Xmid = 0.5*(X(1:end-1,:) + X(2:end,:));

        gradPsi = evaluate_P3_gradient( ...
            P1, P3, psi, Xmid);

        gradMagnitude = hypot( ...
            gradPsi(:,1),gradPsi(:,2));

        if any(~isfinite(gradMagnitude)) || ...
                any(gradMagnitude <= 0)
            error(['Invalid |grad(psi)| on psiN = %.6g. ', ...
                   'Increase the contour resolution.'],psiNk);
        end

        Fk = interp1( ...
            psiNTable,FTable,psiNk,'linear');

        fluxIntegral = sum( ...
            dl./(Xmid(:,1).*gradMagnitude));

        qInterior(k) = ...
            Fk*fluxIntegral/(2*pi);

        minimumGradPsi(k) = min(gradMagnitude);
    end

    % Extrapolate q to psiN = 1.
    [qBoundary, extrapolationSpread] = ...
        extrapolate_boundary_q(psiNLevels,qInterior);

    % Completed q profile including both endpoints.
    qResult = struct();

    qResult.psiN = [0; psiNLevels(:); 1];
    qResult.q = [qAxis; qInterior; qBoundary];

    qResult.qAxis = qAxis;
    qResult.qBoundary = qBoundary;

    qResult.axis = axisInfo;

    qResult.interiorPsiN = psiNLevels(:);
    qResult.interiorQ = qInterior;
    qResult.minimumGradPsi = minimumGradPsi;

    qResult.boundaryExtrapolationSpread = ...
        extrapolationSpread;

    qResult.contourGridSize = [nR,nZ];
end

function [qAxis, info] = evaluate_axis_q_hessian( ...
    P3, psi, axisNode, axisMode, Faxis)
%EVALUATE_AXIS_Q_HESSIAN
% Recover a smooth quadratic representation of psi around the nodal axis,
% refine the axis location, and evaluate
%
%   q_axis = F_axis/(R_axis*sqrt(det(H_psi))).

    axisPointNodal = P3.points(axisNode,:);

    % Elements containing the nodal-axis DOF.
    centralElements = find( ...
        any(P3.elements == axisNode,2));

    if isempty(centralElements)
        error('[Q_AXIS] Failed to locate elements around axisNode.');
    end

    % Construct an expanded patch for Hessian recovery.
    centralNodes = unique( ...
        P3.elements(centralElements,:));

    patchElements = find( ...
        any(ismember(P3.elements,centralNodes),2));

    patchNodes = unique( ...
        P3.elements(patchElements,:));

    point = P3.points(patchNodes,:);
    value = psi(patchNodes);

    displacement = point-axisPointNodal;

    patchScale = max(hypot( ...
        displacement(:,1),displacement(:,2)));

    if ~isfinite(patchScale) || patchScale <= 0
        error('[Q_AXIS] Invalid magnetic-axis patch.');
    end

    % Scaled coordinates improve conditioning.
    x = displacement(:,1)/patchScale;
    y = displacement(:,2)/patchScale;

    % Quadratic model:
    %
    % psi = c1 + c2*x + c3*y
    %       + 0.5*c4*x^2 + c5*x*y + 0.5*c6*y^2
    A = [ ...
        ones(size(x)), ...
        x, y, ...
        0.5*x.^2, x.*y, 0.5*y.^2];

    % Give nearby points greater influence.
    radius = hypot(x,y);
    weight = 1./(1 + (radius/0.5).^4);

    sqrtWeight = sqrt(weight);

    weightedA = A.*sqrtWeight;
    weightedValue = value.*sqrtWeight;

    if rank(weightedA) < 6
        error('[Q_AXIS] Axis patch cannot determine a quadratic Hessian.');
    end

    coefficient = weightedA\weightedValue;

    gradientAtNodalPoint = ...
        [coefficient(2); coefficient(3)]/patchScale;

    Hpsi = [ ...
        coefficient(4), coefficient(5);
        coefficient(5), coefficient(6)]/patchScale^2;

    Hpsi = 0.5*(Hpsi+Hpsi.');

    eigenvalue = eig(Hpsi);

    switch axisMode
        case 'max'
            if any(eigenvalue >= 0)
                error(['[Q_AXIS] Recovered Hessian is not negative ', ...
                       'definite at a maximum-type axis.']);
            end

        case 'min'
            if any(eigenvalue <= 0)
                error(['[Q_AXIS] Recovered Hessian is not positive ', ...
                       'definite at a minimum-type axis.']);
            end

        otherwise
            error('[Q_AXIS] axisMode must be max or min.');
    end

    determinantH = det(Hpsi);

    if ~isfinite(determinantH) || determinantH <= 0
        error('[Q_AXIS] Invalid Hessian determinant.');
    end

    % Stationary point of the recovered quadratic model.
    axisShift = -Hpsi\gradientAtNodalPoint;
    axisPointRecovered = axisPointNodal + axisShift.';

    if norm(axisShift) > 0.5*patchScale
        error(['[Q_AXIS] Recovered axis lies too far from the ', ...
               'nodal-axis patch center.']);
    end

    P1 = triangulation(P3.P1elements,P3.P1points);

    if isnan(pointLocation(P1,axisPointRecovered))
        error('[Q_AXIS] Recovered magnetic axis is outside the mesh.');
    end

    if axisPointRecovered(1) <= 0
        error('[Q_AXIS] Recovered magnetic-axis radius is not positive.');
    end

    qAxis = Faxis/( ...
        axisPointRecovered(1)*sqrt(determinantH));

    % Evaluate recovered quadratic psi at the refined axis.
    u = axisShift(1)/patchScale;
    v = axisShift(2)/patchScale;

    axisBasis = [ ...
        1, u, v, ...
        0.5*u^2, u*v, 0.5*v^2];

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

    info.fitResidual = fitResidual;
    info.patchScale = patchScale;
    info.patchNodeCount = numel(patchNodes);
    info.fitConditionNumber = cond(weightedA);
end

function [qBoundary, spread] = ...
    extrapolate_boundary_q(psiN,q)

    psiN = psiN(:);
    q = q(:);

    if numel(q) < 3
        error('[Q_BOUNDARY] At least three interior q values are required.');
    end

    nFit = min(5,numel(q));
    index = numel(q)-nFit+1:numel(q);

    % Shift the coordinate so the target boundary is x = 0.
    x = psiN(index)-1;
    y = q(index);

    linearCoefficient = polyfit(x,y,1);
    qBoundaryLinear = polyval(linearCoefficient,0);

    if nFit >= 4
        quadraticCoefficient = polyfit(x,y,2);
        qBoundary = polyval(quadraticCoefficient,0);

        spread = abs(qBoundary-qBoundaryLinear);
    else
        qBoundary = qBoundaryLinear;
        spread = NaN;
    end

    if ~isfinite(qBoundary)
        error('[Q_BOUNDARY] Boundary-q extrapolation failed.');
    end

    if isfinite(spread) && ...
            spread > 0.05*max(abs(qBoundary),eps)
        warning('GS:QBoundary:ExtrapolationSensitive', ...
            ['Linear and quadratic q-boundary extrapolations ', ...
             'differ by %.3e.'],spread);
    end
end

function X = select_axis_contour(C, axisPoint, Rlim, Zlim)
%SELECT_AXIS_CONTOUR Select the closed contour enclosing the magnetic axis.

    curves = {};
    column = 1;

    while column <= size(C,2)
        nPoint = round(C(2,column));
        lastColumn = column+nPoint;

        if nPoint >= 3 && lastColumn <= size(C,2)
            curves{end+1} = C(:,column+1:lastColumn).'; %#ok<AGROW>
        end

        column = lastColumn+1;
    end

    domainScale = max([diff(Rlim),diff(Zlim),1]);
    closureTolerance = 1.e-8*domainScale;

    candidates = {};
    areas = [];

    for k = 1:numel(curves)
        curve = curves{k};
        closureGap = norm(curve(end,:)-curve(1,:));

        if closureGap > closureTolerance
            continue
        end

        if closureGap > 10*eps(domainScale)
            curve(end+1,:) = curve(1,:);
        end

        enclosesAxis = inpolygon( ...
            axisPoint(1),axisPoint(2),curve(:,1),curve(:,2));

        if enclosesAxis
            candidates{end+1} = curve; %#ok<AGROW>
            areas(end+1) = polyarea(curve(:,1),curve(:,2)); %#ok<AGROW>
        end
    end

    if isempty(candidates)
        error(['No closed flux surface enclosing the magnetic axis ', ...
               'was found. Increase nR or change psiNLevels.']);
    end

    if numel(candidates) > 1
        warning('GS:QProfile:MultipleAxisContours', ...
            ['Multiple contours enclose the magnetic axis. ', ...
             'The outermost contour will be used.']);
    end

    [~,index] = max(areas);
    X = candidates{index};
end

function gradPsi = evaluate_P3_gradient(P1, P3, psi, Xq)
% Evaluate the elementwise P3 gradient at Cartesian query points.

    elementID = pointLocation(P1, Xq);

    if any(isnan(elementID))
        error('Some flux-surface integration points are outside the mesh.');
    end

    nPoint = size(Xq,1);
    gradPsi = nan(nPoint,2);

    for k = 1:nPoint
        e = elementID(k);

        lambda = cartesianToBarycentric(P1, e, Xq(k,:));
        xi  = lambda(2);
        eta = lambda(3);

        vertexIDs = P3.P1elements(e,:);
        Ve = P3.P1points(vertexIDs,:);

        [~, dNphysical] = eval_P3_element(xi, eta, Ve);

        solutionIDs = P3.elements(e,:);
        gradient = dNphysical*psi(solutionIDs);

        gradPsi(k,:) = gradient.';
    end
end

function result = evaluate_global_parameters( ...
    P3, quad, psi, ...
    psiAxis, dpsi, ...
    psiNTable, pressureTable, ...
    R0, B0, aMinor, Ip)
%EVALUATE_GLOBAL_PARAMETERS
% Compute plasma volume, volume-averaged pressure, beta and li_3.
%
% Axisymmetric volume element:
%
%   dV = 2*pi*R*dR*dZ
%
% Magnetic convention:
%
%   Bp^2 = |grad(psi)|^2/R^2

    mu0 = 4*pi*1.e-7;

    psiNTable = psiNTable(:);
    pressureTable = pressureTable(:);
    psi = psi(:);

    if numel(psiNTable) ~= numel(pressureTable)
        error(['[GLOBAL_PARAMETERS] psiN and pressure profiles ', ...
               'must have equal lengths.']);
    end

    if any(~isfinite(pressureTable))
        error('[GLOBAL_PARAMETERS] Pressure contains NaN or Inf.');
    end

    pressureScale = max([max(abs(pressureTable)),1]);

    if min(pressureTable) < -1.e-10*pressureScale
        error(['[GLOBAL_PARAMETERS] Pressure is negative. ', ...
               'Check whether p0 means pAxis or pBoundary.']);
    end

    if ~isfinite(R0) || R0 <= 0
        error('[GLOBAL_PARAMETERS] R0 must be positive.');
    end

    if ~isfinite(B0) || B0 == 0
        error('[GLOBAL_PARAMETERS] B0 must be finite and nonzero.');
    end

    currentScale = max(1,abs(Ip));

    if ~isfinite(Ip) || abs(Ip) <= eps(currentScale)
        error('[GLOBAL_PARAMETERS] Ip is zero or invalid.');
    end

    L1 = 1-quad.xi-quad.eta;
    L2 = quad.xi;
    L3 = quad.eta;

    plasmaVolume = 0;
    pressureIntegral = 0;
    Bp2VolumeIntegral = 0;

    nElement = size(P3.elements,1);

    for e = 1:nElement
        vertexIDs = P3.P1elements(e,:);
        Ve = P3.P1points(vertexIDs,:);

        R1 = Ve(1,1);
        Z1 = Ve(1,2);

        R2 = Ve(2,1);
        Z2 = Ve(2,2);

        R3 = Ve(3,1);
        Z3 = Ve(3,2);

        J = [ ...
            R2-R1, R3-R1;
            Z2-Z1, Z3-Z1];

        detJ = det(J);

        if detJ <= 0
            error( ...
                '[GLOBAL_PARAMETERS] Invalid element %d.',e);
        end

        Rq = L1*R1+L2*R2+L3*R3;

        if any(Rq <= 0)
            error('[GLOBAL_PARAMETERS] Encountered R <= 0.');
        end

        solutionIDs = P3.elements(e,:);
        psiElement = psi(solutionIDs);

        % Flux at quadrature points
        psiQ = quad.N*psiElement;
        psiNQ = (psiQ-psiAxis)/dpsi;

        % Allow only small FEM overshoots outside [0,1].
        if any(psiNQ < -1.e-3) || any(psiNQ > 1+1.e-3)
            error(['[GLOBAL_PARAMETERS] Quadrature-point psiN ', ...
                   'extends too far outside [0,1].']);
        end

        psiNForProfile = min(max(psiNQ,0),1);

        pressureQ = interp1( ...
            psiNTable,pressureTable, ...
            psiNForProfile,'linear');

        % Physical P3 derivatives
        dNdR = ( ...
              J(2,2)*quad.dNdxi ...
            - J(2,1)*quad.dNdeta)/detJ;

        dNdZ = ( ...
            - J(1,2)*quad.dNdxi ...
            + J(1,1)*quad.dNdeta)/detJ;

        dpsi_dR = dNdR*psiElement;
        dpsi_dZ = dNdZ*psiElement;

        BpSquared = ( ...
            dpsi_dR.^2+dpsi_dZ.^2)./Rq.^2;

        % dA and axisymmetric dV
        dA = quad.w*detJ;
        dV = 2*pi*Rq.*dA;

        plasmaVolume = ...
            plasmaVolume+sum(dV);

        pressureIntegral = ...
            pressureIntegral+sum(pressureQ.*dV);

        Bp2VolumeIntegral = ...
            Bp2VolumeIntegral+sum(BpSquared.*dV);
    end

    pressureVolumeAverage = ...
        pressureIntegral/plasmaVolume;

    betaTor = ...
        2*mu0*pressureVolumeAverage/B0^2;

    betaPol = ...
        4*pressureIntegral/(R0*mu0*Ip^2);

    li3 = ...
        2*Bp2VolumeIntegral/(R0*mu0^2*Ip^2);

    poloidalMagneticEnergy = ...
        Bp2VolumeIntegral/(2*mu0);

    result = struct();

    result.volume = plasmaVolume;
    result.pressureIntegral = pressureIntegral;
    result.pressureVolumeAverage = pressureVolumeAverage;

    result.Bp2VolumeIntegral = Bp2VolumeIntegral;
    result.poloidalMagneticEnergy = poloidalMagneticEnergy;

    result.betaTor = betaTor;
    result.betaPol = betaPol;
    result.li3 = li3;
    result.betaN = (100*betaTor)*aMinor*abs(B0)/(abs(Ip)/1.e6);
end