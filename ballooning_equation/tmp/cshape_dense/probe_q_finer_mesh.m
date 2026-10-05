root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense',sprintf('column_%02d.mat',columnIndex)),'scan');
scan=saved.scan; config=scan.config;
config.Grid.NR=64; config.Grid.NZ=96;
config.QEvaluation.contourNR=400;
config.Continuation.Enabled=false; config.Continuation.MaximumQIterationRestarts=0;
config.PQConstraint.outer.maxIterations=24;
config.Ballooning.NSurfacePoints=4096; config.Ballooning.FourierModes=64;
config.Ballooning.ProjectionIterations=3;
ip=find(~scan.done,1); p0=scan.pressureAmplitude(ip);
warm=scan.predictors{1};
if columnIndex==17
    old=load(fullfile(root,'tmp','cshape_tests','default_validation.mat'),'report');
    warm=old.report.equilibria{2,3};
end
[point,equilibrium,predictor]=evaluate_indented_cshape_scan_point( ...
    config,scan.shapes{1},p0,[],warm);
point.solverOverrides=struct('Grid',config.Grid,'QEvaluation',config.QEvaluation, ...
    'reason','finer native GS mesh at stalled prescribed-q point');
save(fullfile(root,'tmp','cshape_dense',sprintf('probe_q_finer_mesh_%02d.mat',columnIndex)), ...
    'point','equilibrium','predictor','config','ip','-v7.3');
fprintf('FINER MESH I %.3f p0 %.0f beta %.6f complete %d status %s\n', ...
    point.indent,p0,point.betaPercent,point.complete,point.status);
if ~isempty(point.solverDiagnostics), disp(point.solverDiagnostics); end
