root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
data=load(fullfile(root,'tmp','cshape_tests','default_validation.mat'),'report');
config=make_dense_indented_cshape_scan_config;
selection=[1,1;2,1;3,1;3,3];
results=cell(size(selection,1),2);
for j=1:size(selection,1)
    eq=data.report.equilibria{selection(j,1),selection(j,2)};
    for k=1:2
        options=config.Ballooning; options.CellsPerTurn=96*2^(k-1);
        args=namedargs2cell(options); timer=tic;
        results{j,k}=evaluate_equilibrium_ballooning(eq,args{:});
        bal=results{j,k};
        fprintf('PILOT I %.3f beta %.5f cells %d: %s complete %d lambda %.8g (%0.1fs)\n', ...
            data.report.points{selection(j,1),selection(j,2)}.indent, ...
            data.report.points{selection(j,1),selection(j,2)}.betaPercent, ...
            options.CellsPerTurn,bal.status,bal.complete,bal.lambdaMax,toc(timer));
        if ~bal.complete
            first=find(~cellfun(@isempty,bal.failures),1); disp(bal.failures{first});
        end
        save(fullfile(root,'tmp','cshape_dense','pilot_resolution.mat'),'results','selection','config','-v7.3');
    end
end
