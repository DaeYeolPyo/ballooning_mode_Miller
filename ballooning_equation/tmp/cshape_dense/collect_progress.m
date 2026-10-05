root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
[scan,trace,summary]=collect_dense_indented_cshape_beta_scan;
for ii=1:numel(scan.indentation)
    if any(scan.done(:,ii))
        fprintf('COL %02d I %.3f done %d/%d\n',ii,scan.indentation(ii), ...
            sum(scan.done(:,ii)),numel(scan.pressureAmplitude));
    end
end
