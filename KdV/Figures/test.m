%a script to produce plots of the energy change for N = 20

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

N = 20;
epsilon = 0.1;
endtime = 10;
howoften = 10;

%find the exact solution
simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = endtime;   %end of simulation
simulation_params.howoften = howoften;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 256;          %number of positive modes to simulate
simulation_params.initialization = @(x) full_init_KdV(x);  %full simulation

simulation_params.initialization = @(x) complete_init_KdV(x);
simulation_params.order = 4;
simulation_params.N = N;

simulation_params.F_modes
