root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
config=make_dense_indented_cshape_scan_config;
selection=[1,.09;8,.09;1,.30];
results=cell(size(selection,1),1); comparisons=struct([]);
for j=1:size(selection,1)
    saved=load(fullfile(root,'output','indented_cshape_dense', ...
        sprintf('column_%02d.mat',selection(j,1))),'scan');
    target=selection(j,2)*config.BReference^2/(2*(4*pi*1.e-7));
    [~,ip]=min(abs(saved.scan.pressureAmplitude-target));
    eq=saved.scan.equilibria{ip}; p=saved.scan.points{ip};
    assert(~isempty(eq) && p.complete);
    options=config.Ballooning; options.CellsPerTurn=192;
    args=namedargs2cell(options); results{j}=evaluate_equilibrium_ballooning(eq,args{:});
    bal=results{j};
    item=struct('I',p.indent,'betaPercent',p.betaPercent, ...
        'lambda96',p.lambdaMax,'lambda192',bal.lambdaMax,'status96',p.status, ...
        'status192',bal.status,'complete192',bal.complete, ...
        'maximumArrayDifference',max(abs(p.ballooning.lambda-bal.lambda),[],'all'));
    if isempty(comparisons), comparisons=item; else, comparisons(j)=item; end
    fprintf('RESOLUTION I %.3f beta %.6f lambda96 %.8g lambda192 %.8g status %s/%s complete %d\n', ...
        p.indent,p.betaPercent,p.lambdaMax,bal.lambdaMax,p.status,bal.status,bal.complete);
    save(fullfile(root,'tmp','cshape_dense','dense_angular_resolution.mat'), ...
        'comparisons','results','selection','-v7.3');
end
