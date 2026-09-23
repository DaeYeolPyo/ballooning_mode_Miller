function figures = plot_GS_equilibrium(equilibrium)
%PLOT_GS_EQUILIBRIUM Plot fields, profiles, convergence, and residual.
%
%   figures = plot_GS_equilibrium(equilibrium)

    arguments
        equilibrium (1,1) struct
    end

    P3 = equilibrium.mesh;
    result = equilibrium.solverResult;

    % P1 geometry triangulation used to locate plotting points.
    P1 = triangulation(P3.P1elements, P3.P1points);

    Rlim = [min(P3.P1points(:,1)), max(P3.P1points(:,1))];
    Zlim = [min(P3.P1points(:,2)), max(P3.P1points(:,2))];

    % Keep approximately uniform physical sampling in R and Z.
    nR = 220;
    aspectRatio = diff(Zlim)/diff(Rlim);
    nZ = min(500, max(160, round(nR*aspectRatio)));

    Rvec = linspace(Rlim(1), Rlim(2), nR);
    Zvec = linspace(Zlim(1), Zlim(2), nZ);
    [Rplot,Zplot] = meshgrid(Rvec,Zvec);

    queryPoints = [Rplot(:),Zplot(:)];

    % P3-aware field evaluation. Outside points remain NaN.
    psiPlot = eval_P3_sol( ...
        P1, P3, equilibrium.psi, queryPoints);

    psiPlot = reshape(psiPlot,size(Rplot));

    psiNPlot = ...
        (psiPlot-result.psiAxis)/result.dpsi;

    %% Flux fields
    fieldFigure = figure( ...
        'Name','Grad-Shafranov equilibrium', ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','normalized', ...
        'Position',[0.03,0.08,0.94,0.80]);

    layout = tiledlayout(1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    % Dimensional flux
    ax1 = nexttile(layout);

    contourf(ax1,Rplot,Zplot,psiPlot,32, ...
        'LineStyle','none');

    colorbar(ax1);
    title(ax1,'Poloidal flux, \psi');
    decorate_spatial_axes( ...
        ax1,P3,equilibrium.axisPoint,Rlim,Zlim);

    % Normalized flux
    ax2 = nexttile(layout);

    contourf(ax2,Rplot,Zplot,psiNPlot,0:0.05:1, ...
        'LineStyle','none');

    hold(ax2,'on');

    contour(ax2,Rplot,Zplot,psiNPlot,0.1:0.1:0.9, ...
        'k-','LineWidth',0.8);

    clim(ax2,[0,1]);
    colorbar(ax2);
    title(ax2,'Normalized flux, \psi_N');

    decorate_spatial_axes( ...
        ax2,P3,equilibrium.axisPoint,Rlim,Zlim);

    colormap(fieldFigure,turbo);
    title(layout,'Fixed-boundary Grad-Shafranov equilibrium');

    %% Flux profiles
    profile = equilibrium.profile;
    requiredProfileFields = {'psiN','pressure','q','Fpol'};

    for k = 1:numel(requiredProfileFields)
        fieldName = requiredProfileFields{k};

        if ~isfield(profile,fieldName)
            error('GS:Plot:MissingProfile', ...
                'equilibrium.profile.%s is required for plotting.', ...
                fieldName);
        end
    end

    psiNProfile = profile.psiN(:);
    pressure = profile.pressure(:);
    qProfile = profile.q(:);
    Fpol = profile.Fpol(:);

    nProfile = numel(psiNProfile);

    if nProfile < 2 || any(~isfinite(psiNProfile)) || ...
            any(diff(psiNProfile) <= 0)
        error('GS:Plot:InvalidPsiNProfile', ...
            'equilibrium.profile.psiN must be finite and increasing.');
    end

    if numel(pressure) ~= nProfile || numel(qProfile) ~= nProfile || ...
            numel(Fpol) ~= nProfile
        error('GS:Plot:ProfileSize', ...
            'Pressure, q, F, and psiN profiles must have equal lengths.');
    end

    if any(~isfinite(pressure)) || any(~isfinite(qProfile)) || ...
            any(~isfinite(Fpol))
        error('GS:Plot:InvalidProfile', ...
            'Pressure, q, and F profiles must contain finite values.');
    end

    profileFigure = figure( ...
        'Name','Grad-Shafranov profiles', ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','normalized', ...
        'Position',[0.06,0.16,0.88,0.54]);

    profileLayout = tiledlayout(1,3, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    axPressure = nexttile(profileLayout);
    plot(axPressure,psiNProfile,pressure, ...
        'b-','LineWidth',1.5);
    decorate_profile_axes(axPressure,'Pressure [Pa]','Pressure');

    axQ = nexttile(profileLayout);
    plot(axQ,psiNProfile,qProfile, ...
        'r-','LineWidth',1.5);
    decorate_profile_axes(axQ,'q','Safety factor');

    axF = nexttile(profileLayout);
    plot(axF,psiNProfile,Fpol, ...
        'Color',[0.10,0.55,0.20], ...
        'LineWidth',1.5);
    decorate_profile_axes(axF,'F [T m]','Toroidal-field function');

    title(profileLayout,'Flux profiles');

    %% Picard convergence
    convergenceFigure = figure( ...
        'Name','Grad-Shafranov convergence', ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','normalized', ...
        'Position',[0.20,0.16,0.60,0.58]);

    iteration = 1:result.iterations;

    semilogy(iteration,result.updateHistory, ...
        'o-','LineWidth',1.3);

    hold on;

    semilogy(iteration,result.residualHistory, ...
        's-','LineWidth',1.3);

    grid on;
    box on;

    xlabel('Picard iteration');
    ylabel('Relative error');

    legend( ...
        'solution update', ...
        'nonlinear residual', ...
        'Location','best');

    title('Grad-Shafranov Picard convergence');
    set(gca,'FontSize',11);

    %% Final discrete weak residual
    % The entries of K*psi-f are residual functionals evaluated with the
    % P3 basis functions. They are not pointwise strong-form residuals, so
    % plot them at the free-DOF locations rather than interpolating them as
    % a continuous scalar field.
    weakResidualFigure = plot_final_weak_residual( ...
        P3,result,equilibrium.axisPoint,Rlim,Zlim);

    qConstraintFigure = [];

    if isfield(equilibrium,'constraintResult') && ...
            isfield(profile,'qPrescribed') && ...
            ~isempty(profile.qPrescribed)
        qConstraintFigure = plot_pq_q_comparison(equilibrium);
    end

    figures = struct( ...
        'fields',fieldFigure, ...
        'profiles',profileFigure, ...
        'convergence',convergenceFigure, ...
        'weakResidual',weakResidualFigure, ...
        'qConstraint',qConstraintFigure);
end


function decorate_spatial_axes(ax,P3,axisPoint,Rlim,Zlim)
    hold(ax,'on');

    % Draw the polygonal fixed boundary.
    edge = P3.boundaryEdges;
    point = P3.P1points;

    Rboundary = [point(edge(:,1),1), point(edge(:,2),1)].';
    Zboundary = [point(edge(:,1),2), point(edge(:,2),2)].';

    plot(ax,Rboundary,Zboundary, ...
        'k-','LineWidth',1.5);

    plot(ax,axisPoint(1),axisPoint(2), ...
        'rx','MarkerSize',11,'LineWidth',2);

    % A large figure and compact 1x2 layout keep axis-equal plots large.
    axis(ax,'equal');
    xlim(ax,Rlim);
    ylim(ax,Zlim);

    grid(ax,'on');
    box(ax,'on');

    xlabel(ax,'R [m]');
    ylabel(ax,'Z [m]');

    ax.FontSize = 11;
end


function decorate_profile_axes(ax,yLabelText,titleText)
    grid(ax,'on');
    box(ax,'on');

    xlim(ax,[0,1]);

    xlabel(ax,'\psi_N');
    ylabel(ax,yLabelText);
    title(ax,titleText);

    ax.FontSize = 11;
end


function residualFigure = plot_final_weak_residual( ...
    P3,result,axisPoint,Rlim,Zlim)
%PLOT_FINAL_WEAK_RESIDUAL Plot K*psi-f on the unconstrained P3 DOFs.

    requiredFields = {'finalResidual','freeNodes','residualHistory'};

    for iField = 1:numel(requiredFields)
        fieldName = requiredFields{iField};

        if ~isfield(result,fieldName)
            error('GS:Plot:MissingResidualField', ...
                'equilibrium.solverResult.%s is required.',fieldName);
        end
    end

    freeNodes = result.freeNodes(:);
    residual = result.finalResidual(:);
    residualHistory = result.residualHistory(:);
    nDof = size(P3.points,1);

    if isempty(freeNodes) || numel(residual) ~= numel(freeNodes)
        error('GS:Plot:ResidualSize', ...
            ['solverResult.finalResidual must contain one value ', ...
             'per free node.']);
    end

    if any(freeNodes < 1) || any(freeNodes > nDof) || ...
            any(freeNodes ~= round(freeNodes)) || ...
            numel(unique(freeNodes)) ~= numel(freeNodes)
        error('GS:Plot:InvalidFreeNodes', ...
            'solverResult.freeNodes contains invalid node indices.');
    end

    if any(~isfinite(residual)) || isempty(residualHistory) || ...
            any(~isfinite(residualHistory)) || ...
            any(residualHistory < 0)
        error('GS:Plot:InvalidResidual', ...
            'The final residual and residual history must be finite.');
    end

    % Recover the balance scale used by Picard_iteration:
    %
    %   relativeResidual = norm(finalResidual)/balanceScale.
    %
    % This makes the plotted coefficients dimensionless and ensures that
    % their 2-norm equals the reported final relative residual (apart from
    % the exactly-zero case).
    residualNorm = norm(residual);
    finalRelativeResidual = residualHistory(end);

    if residualNorm == 0
        normalizedResidual = zeros(size(residual));
    elseif finalRelativeResidual > 0
        balanceScale = residualNorm/finalRelativeResidual;
        normalizedResidual = residual/balanceScale;
    else
        error('GS:Plot:InconsistentResidual', ...
            ['A nonzero final residual cannot have a zero reported ', ...
             'relative residual.']);
    end

    residualFigure = figure( ...
        'Name','Grad-Shafranov final weak residual', ...
        'NumberTitle','off', ...
        'Color','w', ...
        'Units','normalized', ...
        'Position',[0.08,0.10,0.84,0.70]);

    residualLayout = tiledlayout(residualFigure,1,2, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    % Signed residual coefficients.
    axSigned = nexttile(residualLayout);
    draw_mesh_background(axSigned,P3);

    scatter(axSigned, ...
        P3.points(freeNodes,1), ...
        P3.points(freeNodes,2), ...
        22,normalizedResidual,'filled');

    signedLimit = max(abs(normalizedResidual));

    if signedLimit == 0
        signedLimit = 1;
    end

    clim(axSigned,[-signedLimit,signedLimit]);
    colormap(axSigned,blue_white_red_colormap(256));

    signedColorbar = colorbar(axSigned);
    signedColorbar.Label.String = 'r_i / balance scale';

    title(axSigned,'Signed weak-residual coefficient');
    decorate_spatial_axes(axSigned,P3,axisPoint,Rlim,Zlim);

    % Log-magnitude residual coefficients. The floor prevents log10(0)
    % while remaining far below any practically relevant tolerance.
    axMagnitude = nexttile(residualLayout);
    draw_mesh_background(axMagnitude,P3);

    magnitudeFloor = 10*eps(max(1,max(abs(normalizedResidual))));
    logMagnitude = log10(max(abs(normalizedResidual),magnitudeFloor));

    scatter(axMagnitude, ...
        P3.points(freeNodes,1), ...
        P3.points(freeNodes,2), ...
        22,logMagnitude,'filled');

    colormap(axMagnitude,parula(256));

    magnitudeColorbar = colorbar(axMagnitude);
    magnitudeColorbar.Label.String = ...
        'log_{10}(|r_i| / balance scale)';

    title(axMagnitude,'Weak-residual magnitude');
    decorate_spatial_axes(axMagnitude,P3,axisPoint,Rlim,Zlim);

    title(residualLayout,sprintf( ...
        'Final discrete weak residual, relative 2-norm = %.3e', ...
        norm(normalizedResidual)));
end


function draw_mesh_background(ax,P3)
%DRAW_MESH_BACKGROUND Show the affine P1 element skeleton.

    patch( ...
        'Parent',ax, ...
        'Faces',P3.P1elements, ...
        'Vertices',P3.P1points, ...
        'FaceColor','none', ...
        'EdgeColor',[0.82,0.82,0.82], ...
        'LineWidth',0.35);

    hold(ax,'on');
end


function map = blue_white_red_colormap(nColor)
%BLUE_WHITE_RED_COLORMAP Small dependency-free diverging colormap.

    halfColor = floor(nColor/2);
    lower = [ ...
        linspace(0.10,1.00,halfColor).', ...
        linspace(0.30,1.00,halfColor).', ...
        linspace(0.85,1.00,halfColor).'];

    upperCount = nColor-halfColor;
    upper = [ ...
        linspace(1.00,0.85,upperCount).', ...
        linspace(1.00,0.15,upperCount).', ...
        linspace(1.00,0.10,upperCount).'];

    map = [lower;upper];
end
