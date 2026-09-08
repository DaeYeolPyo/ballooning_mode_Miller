function fig = plot_local_ballooning_coefficients(localBal)
%PLOT_LOCAL_BALLOONING_COEFFICIENTS Plot g, c, f, and K before solving.

narginchk(1,1);
if ~isstruct(localBal) || ~isscalar(localBal) || ...
        ~isfield(localBal,'coefficients')
    error('GS:ballooning:InvalidCoefficientResult', ...
        'localBal must come from build_local_ballooning_coefficients.');
end
coeff = localBal.coefficients;
x = coeff.theta/pi;

fig = figure('Name','Local ballooning coefficients','Color','w');
layout = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
plot(ax,x,coeff.g,'LineWidth',1.5)
grid(ax,'on')
ylabel(ax,'g')
title(ax,'Field-line bending')

ax = nexttile(layout);
plot(ax,x,coeff.c,'LineWidth',1.5)
hold(ax,'on')
yline(ax,0,'k:')
hold(ax,'off')
grid(ax,'on')
ylabel(ax,'c')
title(ax,'Pressure-curvature drive')

ax = nexttile(layout);
plot(ax,x,coeff.f,'LineWidth',1.5)
grid(ax,'on')
xlabel(ax,'\theta/\pi')
ylabel(ax,'f')
title(ax,'Eigenvalue weight')

ax = nexttile(layout);
plot(ax,x,coeff.K,'LineWidth',1.5)
grid(ax,'on')
xlabel(ax,'\theta/\pi')
ylabel(ax,'|\nabla\alpha|^2')
title(ax,'Ballooning wave-vector metric')

title(layout,sprintf( ...
    'Local PEST coefficients at \\psi_N = %.2f (no eigenvalue solve)', ...
    localBal.targetPsiN))
end
