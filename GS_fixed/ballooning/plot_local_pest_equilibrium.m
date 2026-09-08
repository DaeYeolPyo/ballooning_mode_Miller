function fig = plot_local_pest_equilibrium(localEq)
%PLOT_LOCAL_PEST_EQUILIBRIUM Plot the PEST coordinate grid and profiles.

narginchk(1,1);
if ~isstruct(localEq) || ~isscalar(localEq) || ...
        ~all(isfield(localEq,{'map','profiles','bundle'}))
    error('GS:ballooning:InvalidLocalEquilibrium', ...
        'localEq must be returned by build_local_pest_equilibrium.');
end
map = localEq.map;

fig = figure('Name','Local PEST equilibrium','Color','w');
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
hold(ax,'on')
colors = parula(numel(map.psiN));
for k = 1:numel(map.psiN)
    R = map.R(k,[1:end,1]);
    Z = map.Z(k,[1:end,1]);
    plot(ax,R,Z,'-','Color',colors(k,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('\\psi_N = %.2f',map.psiN(k)))
end
stride = max(1,round(numel(map.theta)/16));
indices = 1:stride:numel(map.theta);
for index = indices
    plot(ax,map.R(:,index),map.Z(:,index),'-', ...
        'Color',[0.55,0.55,0.55],'LineWidth',0.8, ...
        'HandleVisibility','off')
end
plot(ax,localEq.bundle.axis.R,localEq.bundle.axis.Z,'kp', ...
    'MarkerFaceColor','y','MarkerSize',10,'DisplayName','magnetic axis')
hold(ax,'off')
axis(ax,'equal','tight')
grid(ax,'on')
xlabel(ax,'R')
ylabel(ax,'Z')
title(ax,'Uniform PEST \theta grid')
legend(ax,'Location','bestoutside')

ax = nexttile(layout);
yyaxis(ax,'left')
plot(ax,localEq.profiles.psiN,localEq.profiles.q,'o-', ...
    'LineWidth',1.5,'DisplayName','q')
hold(ax,'on')
plot(ax,localEq.profiles.psiN,localEq.profiles.sHat,'s--', ...
    'LineWidth',1.5,'DisplayName','$\hat{s}$')
ylabel(ax,'$q,\ \hat{s}$','Interpreter','latex')
yyaxis(ax,'right')
plot(ax,localEq.profiles.psiN,localEq.profiles.alpha,'d-', ...
    'LineWidth',1.5,'DisplayName','$\alpha$')
ylabel(ax,'reference $\alpha$','Interpreter','latex')
xlabel(ax,'\psi_N')
grid(ax,'on')
title(ax,'Local equilibrium profiles')
legend(ax,'Location','best','Interpreter','latex')
end
