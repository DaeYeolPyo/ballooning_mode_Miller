root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
issues=checkcode(fullfile(root,'complete_dense_indented_cshape_column.m'),'-id');
for k=1:numel(issues)
    fprintf('line %d %s: %s\n',issues(k).line,issues(k).id,issues(k).message);
end
