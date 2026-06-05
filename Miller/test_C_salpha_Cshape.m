clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

outC = run_cshape_parameter_salpha_example('C');
