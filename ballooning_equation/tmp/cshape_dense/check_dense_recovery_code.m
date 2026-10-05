root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
files={'complete_dense_indented_cshape_column.m','refine_dense_indented_cshape_geometry.m', ...
    'collect_dense_indented_cshape_beta_scan.m','solve_indented_cshape_equilibrium.m'};
for j=1:numel(files)
    issues=checkcode(fullfile(root,files{j}),'-id');
    fprintf('%s: %d analyzer messages\n',files{j},numel(issues));
    if ~isempty(issues), disp(issues); end
end
[~,~,summary]=collect_dense_indented_cshape_beta_scan;
disp(summary);
