%A script to compute and compare the energy derivative of KdV up to time 10
%against those results as would come from a memory model
%
%Figure 1 in current text

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

simulation_params.name = 'full'; 

%find exact solution
[t_list,u_list] = KdV_solve(simulation_params);

%generate derivative data
simulation_params = full_init(simulation_params);
N_list = 20;
[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);

fig1 = figure;
subplot(3,1,1)
plot(t_list,sum(energy_flow_list(1:20,:))-sum(squeeze(nonlin0_energy_flow(1,:,:))),'k')

title('Exact mass derivative of first $20$ modes','fontsize',16,'interpreter','latex')
xlabel('time','fontsize',16)
ylabel('$\Delta M_F$','fontsize',16,'interpreter','latex')
xlabel('time','fontsize',16,'interpreter','latex')

subplot(3,1,2)
plot(t_list,squeeze(sum(nonlin2_energy_flow(1,:,:))),'k')
title('Mass derivative due to $R_k^2$ of first $20$ modes','fontsize',16,'interpreter','latex')
xlabel('time','fontsize',16)
ylabel('$\Delta M_F^2$','fontsize',16,'interpreter','latex')
xlabel('time','fontsize',16,'interpreter','latex')

subplot(3,1,3)
plot(t_list,squeeze(sum(nonlin4_energy_flow(1,:,:))),'k')
title('Mass derivative due to $R_k^4$ of first $20$ modes','fontsize',16,'interpreter','latex')
xlabel('time','fontsize',16)
ylabel('$\Delta M_F^4$','fontsize',16,'interpreter','latex')
xlabel('time','fontsize',16,'interpreter','latex')

saveas(gcf,'compare_energy','eps')