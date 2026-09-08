function fig = plot_local_salpha_diagram(scan)
%PLOT_LOCAL_SALPHA_DIAGRAM Plot lambda-max and maximizing theta0.

narginchk(1,1);
if ~isstruct(scan) || ~isscalar(scan) || ...
        ~all(isfield(scan,{'sHat','alpha','lambdaMax','bestTheta0'}))
    error('GS:ballooning:InvalidSAlphaScan', ...
        'scan must be returned by scan_local_salpha_diagram.');
end

fig = figure('Name','Local s-alpha stability diagram','Color','w');
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
contourf(ax,scan.alpha,scan.sHat,scan.lambdaMax,24, ...
    'LineStyle','none')
hold(ax,'on')
contour(ax,scan.alpha,scan.sHat,scan.lambdaMax,[0,0], ...
    'k','LineWidth',2)
plot(ax,scan.model.reference.alpha,scan.model.reference.sHat, ...
    'wp','MarkerFaceColor',[0.85,0.15,0.1], ...
    'MarkerEdgeColor','k','MarkerSize',11)
hold(ax,'off')
colorbar(ax)
colormap(ax,turbo)
grid(ax,'on')
xlabel(ax,'\alpha')
ylabel(ax,'$\hat{s}$','Interpreter','latex')
title(ax,'max_{\theta_0}\lambda  (positive: unstable)')

ax = nexttile(layout);
imagesc(ax,scan.alpha,scan.sHat,scan.bestTheta0/pi)
set(ax,'YDir','normal')
hold(ax,'on')
contour(ax,scan.alpha,scan.sHat,scan.lambdaMax,[0,0], ...
    'k','LineWidth',1.5)
plot(ax,scan.model.reference.alpha,scan.model.reference.sHat, ...
    'wp','MarkerFaceColor',[0.85,0.15,0.1], ...
    'MarkerEdgeColor','k','MarkerSize',11)
hold(ax,'off')
colorbar(ax)
colormap(ax,hsv(max(8,numel(scan.theta0))))
clim(ax,[-1,1])
grid(ax,'on')
xlabel(ax,'\alpha')
ylabel(ax,'$\hat{s}$','Interpreter','latex')
title(ax,'maximizing \theta_0/\pi')

title(layout,sprintf( ...
    'Frozen-geometry local PEST scan at \\psi_N=%.2f', ...
    scan.model.targetPsiN))
end
