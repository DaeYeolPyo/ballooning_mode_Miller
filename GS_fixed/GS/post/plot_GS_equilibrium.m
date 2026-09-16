function figures = plot_GS_equilibrium(equilibrium)
%PLOT_GS_EQUILIBRIUM Plot flux fields, flux profiles, and convergence.
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

    figures = struct( ...
        'fields',fieldFigure, ...
        'profiles',profileFigure, ...
        'convergence',convergenceFigure);
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
