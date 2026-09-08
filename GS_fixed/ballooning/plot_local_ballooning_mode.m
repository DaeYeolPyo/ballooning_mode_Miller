function fig = plot_local_ballooning_mode(solution)
%PLOT_LOCAL_BALLOONING_MODE Plot the leading mode and eigenvalue spectrum.

narginchk(1,1);
if ~isstruct(solution) || ~isscalar(solution) || ...
        ~all(isfield(solution,{'thetaDof','lambda','leadingMode'}))
    error('GS:ballooning:InvalidEigenSolution', ...
        'solution must be returned by solve_local_ballooning_mode.');
end

mode = solution.leadingMode/max(abs(solution.leadingMode));
fig = figure('Name','Local ballooning eigenmode','Color','w');
layout = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax = nexttile(layout);
plot(ax,solution.thetaDof/pi,mode,'LineWidth',1.6)
hold(ax,'on')
yline(ax,0,'k:')
hold(ax,'off')
grid(ax,'on')
xlabel(ax,'\theta/\pi')
ylabel(ax,'X/max|X|')
title(ax,'Leading Dirichlet eigenfunction')

ax = nexttile(layout);
plot(ax,1:numel(solution.lambda),solution.lambda,'o-', ...
    'LineWidth',1.4,'MarkerFaceColor',[0.2,0.45,0.8])
hold(ax,'on')
yline(ax,0,'k:','Marginal stability')
hold(ax,'off')
grid(ax,'on')
xlabel(ax,'mode index')
ylabel(ax,'\lambda')
title(ax,'Largest eigenvalues')

if solution.lambda(1)>0
    shortStatus = 'unstable candidate';
else
    shortStatus = 'no positive mode';
end
title(layout,sprintf( ...
    '\\psi_N=%.2f: \\lambda_1=%.6g (%s)', ...
    solution.targetPsiN,solution.lambda(1),shortStatus))
end
