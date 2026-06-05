clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

outG = run_cshape_parameter_salpha_example('G');
