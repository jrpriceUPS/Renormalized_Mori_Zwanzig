%computes t^2-coefficients with least squares fit

clear all;close all

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

%parameter details for simulation
simulation_params.epsilon = 0.1;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = 1e-10;    %tolerance for identifying instabilities
simulation_params.N = 128;        %number of positive modes to simulate


%full model with no approximations
model.name = 'full';
model.renormalize = 0;     %logical indicating we don't have a renormalization step here


[t_list,u_list] = KdV_solve(simulation_params,model);

save t_list t_list
save u_list u_list

simulation_params = full_init(simulation_params);

save simulation_params simulation_params

N_list = 8:4:32;
save N_list N_list

[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow] = generate_deriv_data_4func2(t_list,u_list,simulation_params,N_list);
coeffs_list = no_k_dependence_coeffs2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,1);