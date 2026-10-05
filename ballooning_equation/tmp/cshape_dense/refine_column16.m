root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(root); setup_indented_cshape_scan;
refine_dense_indented_cshape_geometry(fullfile(root,'output','indented_cshape_dense','column_16.mat'));
