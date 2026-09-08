clc
clear
close all

this_dir = fileparts(mfilename('fullpath'));

% Dedicated negative-triangularity validation.  Do not inherit a GEQDSK
% path from the workspace.  The wider radial regression suppresses the
% contour-grid noise in d[asin(delta)]/dr for this equilibrium family.
geqdsk_file = fullfile(this_dir, '..', 'gfiles', 'geqdsk_NT0.6');
miller_fit_dpsi = 8.e-2;
target_psiN = 0.7;

miller_NT_validation = test_Miller( ...
    geqdsk_file, target_psiN, miller_fit_dpsi);

assert(contains(miller_NT_validation.geqdsk_file, 'geqdsk_NT0.6'), ...
    'The NT validation did not use geqdsk_NT0.6.');
assert(miller_NT_validation.param.delta < 0, ...
    'The fitted NT surface does not have negative triangularity.');