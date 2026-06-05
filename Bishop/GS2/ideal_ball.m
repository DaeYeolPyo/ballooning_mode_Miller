clc
clear
close all

%% Initialize
init_ballstab();

%% Run
run_stability_check();
write_stability_ascii_2d();

%% Finish
finish_ballstab();
finish_theta_grid();