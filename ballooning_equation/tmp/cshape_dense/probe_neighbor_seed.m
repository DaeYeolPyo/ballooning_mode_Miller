root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense',sprintf('column_%02d.mat',columnIndex)),'scan');
neighbor=load(fullfile(root,'output','indented_cshape_dense',sprintf('column_%02d.mat',neighborIndex)),'scan');
scan=saved.scan; config=scan.config;
config.Continuation.Enabled=false; config.Continuation.MaximumQIterationRestarts=0;
config.PQConstraint.outer.maxIterations=12;
config.Ballooning.NSurfacePoints=4096; config.Ballooning.FourierModes=64;
config.Ballooning.ProjectionIterations=3;
ip=find(~scan.done,1); p0=scan.pressureAmplitude(ip);
[point,equilibrium,predictor]=evaluate_indented_cshape_scan_point( ...
    config,scan.shapes{1},p0,[],neighbor.scan.equilibria{ip});
point.solverOverrides=struct('initialEquilibrium',sprintf('neighbor column %d at same pressure',neighborIndex));
save(fullfile(root,'tmp','cshape_dense',sprintf('probe_neighbor_seed_%02d.mat',columnIndex)), ...
    'point','equilibrium','predictor','config','ip','-v7.3');
fprintf('NEIGHBOR SEED I %.3f p0 %.0f beta %.6f complete %d status %s\n', ...
    point.indent,p0,point.betaPercent,point.complete,point.status);
if ~isempty(point.solverDiagnostics), disp(point.solverDiagnostics); end
