root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); addpath(fullfile(root,'verification')); setup_indented_cshape_scan;
report=test_indented_cshape_pipeline(RunIntegration=true);
save(fullfile(root,'tmp','cshape_dense','final_pipeline_regression.mat'),'report');
files={'make_dense_indented_cshape_scan_config.m','scan_indented_cshape_beta.m', ...
    'examples/run_dense_indented_cshape_beta_scan.m','collect_dense_indented_cshape_beta_scan.m', ...
    'plot_dense_indented_cshape_beta.m','verification/test_rhs_profile_reuse.m', ...
    '../GS_fixed/GS/fem/assemble_RHS.m','../GS_fixed/GS/fem/construct_RHS.m'};
for j=1:numel(files)
    issues=checkcode(fullfile(root,files{j}),'-id');
    fprintf('%s: %d analyzer messages\n',files{j},numel(issues));
end
