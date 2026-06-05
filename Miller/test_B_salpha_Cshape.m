clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

outB = run_cshape_parameter_salpha_example('B');
