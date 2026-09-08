function fig = plot_local_surface_bundle(bundle)
%PLOT_LOCAL_SURFACE_BUNDLE Plot extracted periodic local flux surfaces.

narginchk(1,1);
if ~isstruct(bundle) || ~isscalar(bundle) || ...
        ~all(isfield(bundle,{'surfaces','axis','psiN'}))
    error('GS:ballooning:InvalidSurfaceBundle', ...
        'bundle must be returned by extract_local_surface_bundle.');
end

fig = figure('Name','Local flux-surface bundle','Color','w');
ax = axes(fig);
hold(ax,'on')
colors = parula(bundle.nSurface);
for k = 1:bundle.nSurface
    surface = bundle.surfaces(k);
    points = surface.points([1:end,1],:);
    plot(ax,points(:,1),points(:,2),'-','Color',colors(k,:), ...
        'LineWidth',1.6, ...
        'DisplayName',sprintf('\\psi_N = %.3f',surface.psiN))
    plot(ax,surface.R(1),surface.Z(1),'o','Color',colors(k,:), ...
        'MarkerFaceColor',colors(k,:),'HandleVisibility','off')
end
plot(ax,bundle.axis.R,bundle.axis.Z,'kp','MarkerFaceColor','y', ...
    'MarkerSize',10,'DisplayName','magnetic axis')
hold(ax,'off')
axis(ax,'equal','tight')
grid(ax,'on')
xlabel(ax,'R')
ylabel(ax,'Z')
title(ax,'FEM flux surfaces resampled at uniform geometric arclength')
legend(ax,'Location','bestoutside')
end
