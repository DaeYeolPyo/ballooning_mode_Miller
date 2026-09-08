%% Miller D-shape GS equilibrium and local s-alpha validation
% This script deliberately reuses the same GS_fixed ballooning pipeline as
% the cosine C-shape. The returned variables are:
%   equilibrium       Miller fixed-boundary GS result
%   millerSAlpha      verified local s-alpha report

main_miller
addpath(fullfile(fileparts(mfilename('fullpath')),'ballooning'))
millerSAlpha = verify_miller_gs_salpha(equilibrium,true);
