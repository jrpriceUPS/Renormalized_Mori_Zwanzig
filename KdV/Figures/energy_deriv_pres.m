%A script to compute and compare the energy derivative of KdV up to time 10
%against those results as would come from a memory model

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

%parameter details for simulation
simulation_params.epsilon = 0.1;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;      %tolerance for identifying instabilities
simulation_params.N = 128;        %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);

simulation_params.initialization = @(x) full_init_KdV(x); 

%find exact solution
[t_list,u_list] = PDE_solve(simulation_params);

%generate derivative data
simulation_params = full_init_KdV(simulation_params);
N_list = 20;
[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);

fig1 = figure;
subplot(2,2,1)
plot(t_list,sum(energy_flow_list(1:20,:))-sum(squeeze(nonlin0_energy_flow(1,:,:))),'b','linewidth',2)
axis([0,10,-.001,.001])
xlabel('time','fontsize',16)
ylabel('Exact dM/dt','fontsize',16)

subplot(2,2,2)
plot(t_list,squeeze(sum(nonlin1_energy_flow(1,:,:))),'b','linewidth',2)
xlabel('time','fontsize',16)
ylabel('dM/dt due to R_1','fontsize',16)
axis([0,10,-10,10])

subplot(2,2,3)
plot(t_list,squeeze(sum(nonlin2_energy_flow(1,:,:))),'b','linewidth',2)
xlabel('time','fontsize',16)
ylabel('dM/dt due to R_2','fontsize',16)
axis([0,10,-10,10])

subplot(2,2,4)
plot(t_list,squeeze(sum(nonlin3_energy_flow(1,:,:))),'b','linewidth',2)
xlabel('time','fontsize',16)
ylabel('dM/dt due to R_3','fontsize',16)
axis([0,10,-100000,100000])

saveas(gcf,'energy_deriv','png')