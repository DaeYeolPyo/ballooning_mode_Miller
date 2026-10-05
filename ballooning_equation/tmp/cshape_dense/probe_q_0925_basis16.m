root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense','column_14.mat'),'scan');
scan=saved.scan; config=scan.config; config.PQConstraint.FFprime.nBasis=16;
config.PQConstraint.FFprime.finiteDifferenceAbsoluteStep=.005;
config.PQConstraint.FFprime.finiteDifferenceRelativeStep=.002;
config.PQConstraint.outer.maxIterations=12;
config.Continuation.MaxAttempts=4; config.Continuation.MaximumQIterationRestarts=0;
ip=find(~scan.done,1); p0=scan.pressureAmplitude(ip);
[point,equilibrium,predictor]=evaluate_indented_cshape_scan_point( ...
    config,scan.shapes{1},p0,[],scan.predictors{1});
point.solverOverrides=struct('nBasis',16,'finiteDifferenceAbsoluteStep',.005, ...
    'finiteDifferenceRelativeStep',.002,'reason','q constraint recovery with more profile freedom');
save(fullfile(root,'tmp','cshape_dense','probe_q_0925_basis16.mat'), ...
    'point','equilibrium','predictor','config','ip','-v7.3');
fprintf('BASIS16 RECOVERY p0 %.0f beta %.6f complete %d status %s\n', ...
    p0,point.betaPercent,point.complete,point.status);
if ~isempty(point.solverDiagnostics), disp(point.solverDiagnostics); end
