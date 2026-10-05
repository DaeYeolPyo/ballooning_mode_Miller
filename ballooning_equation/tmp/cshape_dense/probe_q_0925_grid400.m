root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense','column_14.mat'),'scan');
scan=saved.scan; config=scan.config; config.QEvaluation.contourNR=400;
config.PQConstraint.FFprime.finiteDifferenceAbsoluteStep=.005;
config.PQConstraint.FFprime.finiteDifferenceRelativeStep=.002;
config.Continuation.Enabled=false; config.Continuation.MaximumQIterationRestarts=0;
config.PQConstraint.outer.maxIterations=12;
ip=find(~scan.done,1); p0=scan.pressureAmplitude(ip);
[point,equilibrium,predictor]=evaluate_indented_cshape_scan_point( ...
    config,scan.shapes{1},p0,[],scan.predictors{1});
point.solverOverrides=struct('qContourNR',400,'finiteDifferenceAbsoluteStep',.005, ...
    'finiteDifferenceRelativeStep',.002,'reason','q-contour resolution recovery');
save(fullfile(root,'tmp','cshape_dense','probe_q_0925_grid400.mat'), ...
    'point','equilibrium','predictor','config','ip','-v7.3');
fprintf('GRID400 RECOVERY p0 %.0f beta %.6f complete %d status %s\n', ...
    p0,point.betaPercent,point.complete,point.status);
if ~isempty(point.solverDiagnostics), disp(point.solverDiagnostics); end
