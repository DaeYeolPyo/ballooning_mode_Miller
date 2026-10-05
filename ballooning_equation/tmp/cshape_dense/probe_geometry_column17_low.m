root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense','column_17.mat'),'scan');
equilibrium=saved.scan.equilibria{5}; point=saved.scan.points{5};
options=saved.scan.config.Ballooning;
options.NSurfacePoints=4096; options.FourierModes=64; options.ProjectionIterations=3;
args=namedargs2cell(options); result=evaluate_equilibrium_ballooning(equilibrium,args{:});
fprintf('LOW GEOMETRY I=1 beta %.6f: complete %d status %s lambda %.9g\n', ...
    point.betaPercent,result.complete,result.status,result.lambdaMax);
if ~result.complete
    first=find(~cellfun(@isempty,result.failures),1); disp(result.failures{first});
end
save(fullfile(root,'tmp','cshape_dense','probe_geometry_column17_low.mat'), ...
    'result','equilibrium','point','options','-v7.3');
