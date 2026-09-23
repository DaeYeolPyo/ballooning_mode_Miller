function comparison = plot_pq_route_comparison( ...
    equilibriumMode1, equilibriumMode3, options)
%PLOT_PQ_ROUTE_COMPARISON Compare direct and inverse GS solution routes.
%
%   comparison = plot_pq_route_comparison(eqMode1,eqMode3)
%   comparison = plot_pq_route_comparison(eqMode1,eqMode3,options)
%
% The figure contains:
%   1. overlaid normalized-flux surfaces
%   2. spatial difference in normalized flux
%   3. q profiles
%   4. F profiles
%   5. FF' profiles
%   6. normalized profile-difference magnitudes
%
% Options:
%   visible       'on' or 'off'. Default: 'on'
%   gridSize      Cartesian comparison-grid size. Default: 240
%   contourLevels Normalized-flux contour levels. Default: 0.1:0.1:0.9
%   savePath      Optional image path passed to exportgraphics.

    narginchk(2, 3);

    if nargin < 3 || isempty(options)
        options = struct();
    end

    validate_equilibrium(equilibriumMode1,'equilibriumMode1');
    validate_equilibrium(equilibriumMode3,'equilibriumMode3');

    visible = get_option(options,'visible','on');
    gridSize = get_option(options,'gridSize',240);
    contourLevels = get_option( ...
        options,'contourLevels',(0.1:0.1:0.9).');
    savePath = get_option(options,'savePath','');

    if isstring(visible)
        visible = char(visible);
    end

    if ~ismember(visible,{'on','off'})
        error('GS:PQComparison:Visible', ...
            'options.visible must be ''on'' or ''off''.');
    end

    validateattributes(gridSize,{'numeric'}, ...
        {'real','finite','scalar','integer','>=',80});

    contourLevels = contourLevels(:);

    if isempty(contourLevels) || any(~isfinite(contourLevels)) || ...
            any(contourLevels <= 0) || any(contourLevels >= 1)
        error('GS:PQComparison:ContourLevels', ...
            'contourLevels must be finite and lie in (0,1).');
    end

    mesh1 = equilibriumMode1.mesh;
    mesh3 = equilibriumMode3.mesh;

    if ~isequal(mesh1.P1elements,mesh3.P1elements) || ...
            ~isequal(size(mesh1.P1points),size(mesh3.P1points)) || ...
            norm(mesh1.P1points-mesh3.P1points,inf) > 1.0e-12
        error('GS:PQComparison:MeshMismatch', ...
            'Mode-1 and mode-3 equilibria must use the same mesh.');
    end

    vertexIDs = mesh1.vertices(:);
    point = mesh1.P1points;
    psiN1Vertex = equilibriumMode1.psiN(vertexIDs);
    psiN3Vertex = equilibriumMode3.psiN(vertexIDs);

    boundaryPoint = mesh1.P1points( ...
        unique(mesh1.boundaryEdges(:)),:);

    Rvector = linspace(min(point(:,1)),max(point(:,1)),gridSize);
    Zspan = max(point(:,2))-min(point(:,2));
    Rspan = max(point(:,1))-min(point(:,1));
    nZ = max(80,round(gridSize*Zspan/max(Rspan,eps)));
    Zvector = linspace(min(point(:,2)),max(point(:,2)),nZ);
    [Rgrid,Zgrid] = meshgrid(Rvector,Zvector);

    interpolant1 = scatteredInterpolant( ...
        point(:,1),point(:,2),psiN1Vertex,'linear','none');
    interpolant3 = scatteredInterpolant( ...
        point(:,1),point(:,2),psiN3Vertex,'linear','none');

    psiN1Grid = interpolant1(Rgrid,Zgrid);
    psiN3Grid = interpolant3(Rgrid,Zgrid);

    boundaryLoop = order_boundary_loop(mesh1);
    inside = inpolygon( ...
        Rgrid,Zgrid,boundaryLoop(:,1),boundaryLoop(:,2));

    psiN1Grid(~inside) = NaN;
    psiN3Grid(~inside) = NaN;
    deltaPsiNGrid = psiN3Grid-psiN1Grid;

    [profilePsiN,profile1,profile3] = ...
        interpolate_profiles(equilibriumMode1,equilibriumMode3);

    difference = struct();
    difference.q = normalized_difference(profile3.q,profile1.q);
    difference.F = normalized_difference(profile3.F,profile1.F);
    difference.FFprime = normalized_difference( ...
        profile3.FFprime,profile1.FFprime);

    metrics = struct();
    metrics.psiNRootMeanSquare = sqrt(mean( ...
        (equilibriumMode3.psiN-equilibriumMode1.psiN).^2));
    metrics.psiNMaximum = norm( ...
        equilibriumMode3.psiN-equilibriumMode1.psiN,inf);
    metrics.qRelativeRootMeanSquare = relative_rms( ...
        profile3.q,profile1.q);
    metrics.qRelativeMaximum = relative_maximum( ...
        profile3.q,profile1.q);
    metrics.FRelativeL2 = relative_l2(profile3.F,profile1.F);
    metrics.FFprimeRelativeL2 = relative_l2( ...
        profile3.FFprime,profile1.FFprime);

    figureHandle = figure( ...
        'Color','w', ...
        'Visible',visible, ...
        'Name','Mode 1 and mode 3 GS route comparison');

    layout = tiledlayout(figureHandle,2,3, ...
        'TileSpacing','compact','Padding','compact');
    axisHandles = gobjects(6,1);

    %==============================================================
    % Flux-surface overlay
    %==============================================================
    axisHandles(1) = nexttile(layout,1);
    axisHandle = axisHandles(1);
    [~,mode1Contour] = contour( ...
        axisHandle,Rgrid,Zgrid,psiN1Grid, ...
        contourLevels,'b-','LineWidth',1.25);
    hold(axisHandle,'on');
    [~,mode3Contour] = contour( ...
        axisHandle,Rgrid,Zgrid,psiN3Grid, ...
        contourLevels,'r--','LineWidth',1.25);
    boundaryHandle = plot( ...
        axisHandle,boundaryLoop(:,1),boundaryLoop(:,2),'k-');
    mode1AxisHandle = plot( ...
        axisHandle,equilibriumMode1.axisPoint(1), ...
        equilibriumMode1.axisPoint(2),'bo','MarkerFaceColor','b');
    mode3AxisHandle = plot( ...
        axisHandle,equilibriumMode3.axisPoint(1), ...
        equilibriumMode3.axisPoint(2),'rx','LineWidth',1.5);
    axis(axisHandle,'equal');
    axis(axisHandle,'tight');
    grid(axisHandle,'on');
    xlabel(axisHandle,'R [m]');
    ylabel(axisHandle,'Z [m]');
    title(axisHandle,'Normalized-flux surfaces');
    legend(axisHandle, ...
        [mode1Contour,mode3Contour,boundaryHandle, ...
         mode1AxisHandle,mode3AxisHandle], ...
        {'mode 1','mode 3','boundary','axis mode 1','axis mode 3'}, ...
        'Location','best');

    %==============================================================
    % Spatial normalized-flux difference
    %==============================================================
    axisHandles(2) = nexttile(layout,2);
    axisHandle = axisHandles(2);
    contourf(axisHandle,Rgrid,Zgrid,deltaPsiNGrid,24, ...
        'LineStyle','none');
    hold(axisHandle,'on');
    plot(axisHandle,boundaryLoop(:,1),boundaryLoop(:,2),'k-');
    axis(axisHandle,'equal');
    axis(axisHandle,'tight');
    xlabel(axisHandle,'R [m]');
    ylabel(axisHandle,'Z [m]');
    title(axisHandle,'\Delta\psi_N = mode 3 - mode 1');
    colorbar(axisHandle);
    colormap(axisHandle,bluewhitered_colormap(257));

    %==============================================================
    % q, F, and FF' profiles
    %==============================================================
    axisHandles(3) = nexttile(layout,3);
    plot_profile(axisHandles(3),profilePsiN, ...
        profile1.q,profile3.q,'q','Safety factor');

    axisHandles(4) = nexttile(layout,4);
    plot_profile(axisHandles(4),profilePsiN, ...
        profile1.F,profile3.F,'F [T m]','Toroidal-field function');

    axisHandles(5) = nexttile(layout,5);
    plot_profile(axisHandles(5),profilePsiN, ...
        profile1.FFprime,profile3.FFprime, ...
        'F dF/d\psi_N [(T m)^2]','GS source profile');

    %==============================================================
    % Normalized difference magnitudes
    %==============================================================
    axisHandles(6) = nexttile(layout,6);
    axisHandle = axisHandles(6);
    semilogy(axisHandle,profilePsiN, ...
        max(abs(difference.q),eps),'LineWidth',1.3);
    hold(axisHandle,'on');
    semilogy(axisHandle,profilePsiN, ...
        max(abs(difference.F),eps),'LineWidth',1.3);
    semilogy(axisHandle,profilePsiN, ...
        max(abs(difference.FFprime),eps),'LineWidth',1.3);
    grid(axisHandle,'on');
    box(axisHandle,'on');
    xlim(axisHandle,[0,1]);
    xlabel(axisHandle,'\psi_N');
    ylabel(axisHandle,'|mode 3 - mode 1| / max|mode 1|');
    title(axisHandle,'Normalized profile differences');
    legend(axisHandle,{'q','F','FF'''},'Location','best');

    title(layout,sprintf([ ...
        'Mode 1 vs mode 3: RMS(\\Delta\\psi_N)=%.3e, ', ...
        'RMS(\\Delta q/q)=%.3e, L2(\\Delta FF'')=%.3e'], ...
        metrics.psiNRootMeanSquare, ...
        metrics.qRelativeRootMeanSquare, ...
        metrics.FFprimeRelativeL2));

    if ~isempty(savePath)
        exportgraphics(figureHandle,savePath,'Resolution',180);
    end

    comparison = struct();
    comparison.figure = figureHandle;
    comparison.layout = layout;
    comparison.axes = axisHandles;
    comparison.metrics = metrics;
    comparison.profilePsiN = profilePsiN;
    comparison.profileMode1 = profile1;
    comparison.profileMode3 = profile3;
    comparison.normalizedDifference = difference;
    comparison.grid = struct( ...
        'R',Rgrid,'Z',Zgrid,'deltaPsiN',deltaPsiNGrid);
    comparison.boundaryPointCloud = boundaryPoint;
end


function plot_profile(axisHandle,psiN,value1,value3,yLabelText,titleText)
%PLOT_PROFILE Plot a pair of route profiles with common styling.

    plot(axisHandle,psiN,value1,'b-','LineWidth',1.5);
    hold(axisHandle,'on');
    plot(axisHandle,psiN,value3,'r--','LineWidth',1.5);
    grid(axisHandle,'on');
    box(axisHandle,'on');
    xlim(axisHandle,[0,1]);
    xlabel(axisHandle,'\psi_N');
    ylabel(axisHandle,yLabelText);
    title(axisHandle,titleText);
    legend(axisHandle,{'mode 1','mode 3'},'Location','best');
end


function boundaryLoop = order_boundary_loop(mesh)
%ORDER_BOUNDARY_LOOP Follow the P1 boundary edges into a closed polygon.

    edges = mesh.boundaryEdges;
    nEdge = size(edges,1);
    used = false(nEdge,1);
    vertexOrder = zeros(nEdge+1,1);
    vertexOrder(1:2) = edges(1,:).';
    used(1) = true;

    for k = 2:nEdge
        currentVertex = vertexOrder(k);
        candidate = find(~used & ...
            (edges(:,1) == currentVertex | ...
             edges(:,2) == currentVertex),1,'first');

        if isempty(candidate)
            error('GS:PQComparison:BoundaryTopology', ...
                'Failed to order the P1 boundary edges.');
        end

        used(candidate) = true;

        if edges(candidate,1) == currentVertex
            vertexOrder(k+1) = edges(candidate,2);
        else
            vertexOrder(k+1) = edges(candidate,1);
        end
    end

    if vertexOrder(end) ~= vertexOrder(1)
        error('GS:PQComparison:BoundaryTopology', ...
            'The P1 boundary edge chain is not closed.');
    end

    boundaryLoop = mesh.P1points(vertexOrder,:);
end


function [psiN,profile1,profile3] = ...
    interpolate_profiles(equilibrium1,equilibrium3)
%INTERPOLATE_PROFILES Put both routes on the mode-1 profile grid.

    psiN = equilibrium1.profile.psiN(:);
    psiN3 = equilibrium3.profile.psiN(:);

    profile1 = struct();
    profile3 = struct();

    profile1.q = equilibrium1.profile.q(:);
    profile1.F = equilibrium1.profile.Fpol(:);
    profile1.G = equilibrium1.profile.Gpol(:);
    profile1.FFprime = equilibrium1.profile.FdF_dpsiN(:);

    profile3.q = interp1( ...
        psiN3,equilibrium3.profile.q(:),psiN,'pchip');
    profile3.F = interp1( ...
        psiN3,equilibrium3.profile.Fpol(:),psiN,'pchip');
    profile3.G = interp1( ...
        psiN3,equilibrium3.profile.Gpol(:),psiN,'pchip');
    profile3.FFprime = interp1( ...
        psiN3,equilibrium3.profile.FdF_dpsiN(:),psiN,'pchip');
end


function difference = normalized_difference(candidate,reference)
%NORMALIZED_DIFFERENCE Difference normalized by a global reference scale.

    referenceScale = max(norm(reference,inf),eps);
    difference = (candidate-reference)/referenceScale;
end


function value = relative_rms(candidate,reference)
%RELATIVE_RMS Pointwise relative RMS with a scale-aware denominator.

    referenceScale = max(norm(reference,inf),1);
    denominator = max(abs(reference),1.0e-12*referenceScale);
    value = sqrt(mean(((candidate-reference)./denominator).^2));
end


function value = relative_maximum(candidate,reference)
%RELATIVE_MAXIMUM Maximum pointwise relative difference.

    referenceScale = max(norm(reference,inf),1);
    denominator = max(abs(reference),1.0e-12*referenceScale);
    value = max(abs((candidate-reference)./denominator));
end


function value = relative_l2(candidate,reference)
%RELATIVE_L2 Euclidean relative error.

    value = norm(candidate-reference)/max(norm(reference),eps);
end


function colorMap = bluewhitered_colormap(nColor)
%BLUEWHITERED_COLORMAP Symmetric blue-white-red difference colors.

    nLower = ceil(nColor/2);
    nUpper = floor(nColor/2);
    lower = [ ...
        linspace(0.12,1,nLower).', ...
        linspace(0.32,1,nLower).', ...
        ones(nLower,1)];
    upper = [ ...
        ones(nUpper,1), ...
        linspace(1,0.20,nUpper).', ...
        linspace(1,0.16,nUpper).'];
    colorMap = [lower;upper];
end


function validate_equilibrium(equilibrium,valueName)
%VALIDATE_EQUILIBRIUM Validate fields needed by the comparison figure.

    requiredFields = {'mesh','psiN','axisPoint','profile'};

    if ~isstruct(equilibrium) || ~isscalar(equilibrium)
        error('GS:PQComparison:Equilibrium', ...
            '%s must be a scalar structure.',valueName);
    end

    for k = 1:numel(requiredFields)
        if ~isfield(equilibrium,requiredFields{k})
            error('GS:PQComparison:MissingField', ...
                '%s.%s is required.',valueName,requiredFields{k});
        end
    end

    profileFields = {'psiN','q','Fpol','Gpol','FdF_dpsiN'};

    for k = 1:numel(profileFields)
        if ~isfield(equilibrium.profile,profileFields{k})
            error('GS:PQComparison:MissingProfileField', ...
                '%s.profile.%s is required.',valueName,profileFields{k});
        end
    end
end


function value = get_option(options,fieldName,defaultValue)
%GET_OPTION Read an optional structure field.

    if isfield(options,fieldName) && ~isempty(options.(fieldName))
        value = options.(fieldName);
    else
        value = defaultValue;
    end
end
