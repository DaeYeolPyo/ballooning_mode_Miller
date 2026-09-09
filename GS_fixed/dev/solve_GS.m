function equilibrium = solve_GS(input)
    totalTimer = tic;
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
    R0 = input.dim.R0;
    
    % Profile inputs
    % To be developed
    % profile_mode = 1: p' and FF'
    %              = 2: p' and <J_par>
    %              = 3: p' and q
    profile_mode = input.profile.mode;
    pprime       = input.profile.pprime;
    FFprime      = input.profile.FFprime;
    p0           = input.profile.p0;

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
    verbose = input.output.verbose;
    checkTime = input.output.checkTime;

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

    plasmaBoundary = [R_bnd; Z_bnd];
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
        error(['[SOLVE_GS] The constrianed triangulation' ...
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
    %% D.1 Magnetic field and toroidal current evaluation
    psiNProfile = linspace(0, 1, nProfile);
    pressure  = cumtrapz(psiNProfile, pprime) + p0;
    Fpol     = sqrt(cumtrapz(psiNProfile, 2*FFprime));
    Baxis    = Fpol(1)/axisPoint(1);

    %% E.1 Equilibrium result wrap-up
    equilibrium = struct();

    equilibrium.mesh = P3;

    equilibrium.psi = psi;
    equilibrium.psiN = psiN;

    equilibrium.axisNode = axisNode;

    equilibrium.axisPoint = axisPoint;

    equilibrium.profile = struct();
    
    equilibrium.profile.psiN = psiNProfile;
    
    equilibrium.profile.pressure = pressure;
    
    equilibrium.profile.dp_dpsiN = pprime;
    
    equilibrium.profile.Fpol = Fpol;
    
    equilibrium.profile.FdF_dpsiN = FFprime;
      
    equilibrium.dimensionless = struct();

    equilibrium.dimensionless.Raxis = axisPoint(1);
    equilibrium.dimensionless.Zaxis = axisPoint(2);
    equilibrium.dimensionless.simag = psiAxis;
    equilibrium.dimensionless.sibry = psiBoundary;
    equilibrium.dimensionless.Baxis = Baxis;
    %equilibrium.dimensionless.Ip    = Ip;
    
    equilibrium.solverResult = PicardResult;
    %% E.1 Clipping for profile visualization
%     psiNForPlot = min(max(psiN, 0), 1);
%     pressureNodal = interp1(psiNProfile, pressureProfile, ...
%         psiNForPlot, 'linear');
%     F2Nodal = interp1(psiNprofile, F2Profile, ...
%         psiNForPlot, 'linear');
%     FNodal = sqrt(F2Nodal);
%     BphiNodal = Fnodal./P3.points(:, 1);

end