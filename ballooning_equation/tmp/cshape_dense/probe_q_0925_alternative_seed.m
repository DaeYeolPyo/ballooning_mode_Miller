root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense','column_14.mat'),'scan');
old=load(fullfile(root,'tmp','cshape_tests','default_validation.mat'),'report');
scan=saved.scan; config=scan.config;
config.Continuation.Enabled=false; config.Continuation.MaximumQIterationRestarts=0;
config.PQConstraint.outer.maxIterations=12;
ip=find(~scan.done,1); p0=scan.pressureAmplitude(ip);
warm=old.report.equilibria{2,3};
[point,equilibrium,predictor]=evaluate_indented_cshape_scan_point( ...
    config,scan.shapes{1},p0,[],warm);
point.solverOverrides=struct('initialEquilibrium','validated I=1 beta0=0.1 equilibrium', ...
    'reason','alternative converged initial profile at a stalled pressure');
save(fullfile(root,'tmp','cshape_dense','probe_q_0925_alternative_seed.mat'), ...
    'point','equilibrium','predictor','config','ip','-v7.3');
fprintf('ALTERNATIVE SEED p0 %.0f beta %.6f complete %d status %s\n', ...
    p0,point.betaPercent,point.complete,point.status);
if ~isempty(point.solverDiagnostics), disp(point.solverDiagnostics); end
