root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense','column_03.mat'),'scan');
warm=saved.scan.equilibria{find(saved.scan.done,1,'last')};
config=make_dense_indented_cshape_scan_config;
config.PQConstraint.FFprime.nBasis=12;
shape=saved.scan.shapes{1}; p0=.05*config.BReference^2/(2*(4*pi*1.e-7));
timer=tic;
[equilibrium,diagnostics]=solve_indented_cshape_equilibrium(config,shape,p0,[],warm);
wallTime=toc(timer);
save(fullfile(root,'tmp','cshape_dense','probe_basis12.mat'),'equilibrium','diagnostics','config','wallTime','-v7.3');
fprintf('PROBE basis12 I %.3f p0 %.0f success %d time %.1fs\n',shape.indent,p0,~isempty(equilibrium),wallTime);
if ~isempty(equilibrium)
    fprintf('q max %.8g RMS %.8g iterations %d\n', ...
        equilibrium.constraintResult.finalMaxAbsRelativeQ, ...
        equilibrium.constraintResult.finalRmsLogQ,equilibrium.constraintResult.outerIterations);
end
