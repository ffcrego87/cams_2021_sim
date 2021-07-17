%% Quantized distributed estimation

clc;close all;
clearvars

addpath('./util')

%% Problem definition
problem

%% Setup
setup_script

%% Simulation
if Nb_test
    recursive_sim
else
    simulation
end

%% Plot
if Nb_test
    plot_nb_test
else
    plot_sim
end