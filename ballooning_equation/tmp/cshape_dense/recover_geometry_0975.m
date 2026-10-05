root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
saved=load(fullfile(root,'output','indented_cshape_dense','column_16.mat'),'scan');
equilibrium=saved.scan.equilibria{1}; point=saved.scan.points{1};
config=saved.scan.config; results=cell(2,1);
for j=1:2
    options=config.Ballooning;
    options.NSurfacePoints=4096; options.FourierModes=64;
    options.ProjectionIterations=3;
    options.AdapterNR=513+287*(j-1);
    args=namedargs2cell(options); results{j}=evaluate_equilibrium_ballooning(equilibrium,args{:});
    bal=results{j};
    fprintf('GEOMETRY RECOVERY I %.3f beta %.6f adapter %d modes %d: complete %d status %s lambda %.8g\n', ...
        point.indent,point.betaPercent,options.AdapterNR,options.FourierModes, ...
        bal.complete,bal.status,bal.lambdaMax);
    if ~bal.complete
        first=find(~cellfun(@isempty,bal.failures),1); disp(bal.failures{first});
    end
    save(fullfile(root,'tmp','cshape_dense','recover_geometry_0975.mat'), ...
        'results','equilibrium','point','config','-v7.3');
    if bal.complete, break; end
end
