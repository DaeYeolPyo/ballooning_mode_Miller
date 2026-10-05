root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
file=fullfile(root,'output','indented_cshape_dense','column_14.mat');
saved=load(file,'scan'); scan=saved.scan; config=scan.config;
config.PQConstraint.FFprime.finiteDifferenceAbsoluteStep=.005;
config.PQConstraint.FFprime.finiteDifferenceRelativeStep=.002;
ip=find(~scan.done,1); p0=scan.pressureAmplitude(ip);
shape=scan.shapes{1}; warm=scan.predictors{1};
[point,equilibrium,predictor]=evaluate_indented_cshape_scan_point(config,shape,p0,[],warm);
point.solverOverrides=struct('finiteDifferenceAbsoluteStep',.005, ...
    'finiteDifferenceRelativeStep',.002,'reason','smaller q-response perturbations at a stalled pressure');
save(fullfile(root,'tmp','cshape_dense','recover_q_0925.mat'), ...
    'point','equilibrium','predictor','config','ip','-v7.3');
fprintf('RECOVERY I %.3f p0 %.0f beta %.6f complete %d status %s\n', ...
    shape.indent,p0,point.betaPercent,point.complete,point.status);
if point.complete && ~isempty(equilibrium)
    scan.points{ip}=point; scan.equilibria{ip}=equilibrium; scan.predictors{1}=predictor;
    scan.betaPercent(ip)=point.betaPercent; scan.lambdaMax(ip)=point.lambdaMax;
    scan.status(ip)=point.status; scan.done(ip)=true;
    temporary=[tempname(fileparts(file)),'.mat'];
    save(temporary,'scan','-v7.3'); movefile(temporary,file,'f');
end
