%clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

simulation_params.R = inf;          %Reynolds number
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;      %tolerance for identifying instabilities
simulation_params.N = 32;        %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);  %initial condition



%full model with no approximations
simulation_params.initialization = @(x) full_init_BKdV(x);

%non-dimensionalization

U = sqrt(1/(2*pi)*integral(@(x) simulation_params.initial_condition(x).^2,0,2*pi));
simulation_params.epsilon = sqrt(U);
simulation_params = nondim_BKdV(simulation_params);

[t_list,u_list] = PDE_solve(simulation_params);
t_list = t_list * simulation_params.t_scaling;

figure
plot(t_list,get_energy(u_list,8))

pause
%single_plot(t_list(1:1000:end),u_list(:,1:1000:end),32)
single_plot(t_list,u_list,32)