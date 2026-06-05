clc
clear
close all

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);

outH = run_cshape_parameter_salpha_example('H');
