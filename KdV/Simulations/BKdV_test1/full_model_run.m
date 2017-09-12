clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

simulation_params.R = 10;          %Reynolds number
simulation_params.dt = 1e-6;      %timestep
simulation_params.endtime = .1;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;      %tolerance for identifying instabilities
simulation_params.N = 32;        %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);  %initial condition


%full model with no approximations
simulation_params.initialization = @(x) full_init_BKdV(x);

[t_list,u_list] = PDE_solve(simulation_params);

figure
plot(t_list,get_energy(u_list,8))

pause
single_plot(t_list,u_list,32)