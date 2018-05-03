%shows how exact energy derivative compares to ROM energy derivatives R1 through R3
addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

%parameter details for simulation
simulation_params.epsilon = 0.1;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 10;  %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;      %tolerance for identifying instabilities
simulation_params.N = 128;        %number of positive modes to simulate

simulation_params.initialization = @(x) full_init_KdV(x);   

%find exact solution
[t_list,u_list] = PDE_solve(simulation_params);

simulation_params = full_init_KdV(simulation_params);

[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,20);


energy_deriv = figure;

subplot(2,2,1)
set(gca,'FontSize',16)
plot(t_list,sum(energy_flow_list(1:20,:)),'linewidth',2)
xlabel('time')
ylabel('Exact dM/dt')
axis([0,10,-1e-3,1e-3])


subplot(2,2,2)
set(gca,'FontSize',16)
plot(t_list,sum(squeeze(nonlin1_energy_flow(1,1:20,:))),'linewidth',2)
xlabel('time')
ylabel('dM/dt due to R_1')
axis([0,10,-10,10])


subplot(2,2,3)
set(gca,'FontSize',16)
plot(t_list,sum(squeeze(nonlin2_energy_flow(1,1:20,:))),'linewidth',2)
xlabel('time')
ylabel('dM/dt due to R_2')
axis([0,10,-10,10])


subplot(2,2,4)
set(gca,'FontSize',16)
plot(t_list,sum(squeeze(nonlin3_energy_flow(1,1:20,:))),'linewidth',2)
xlabel('time')
ylabel('dM/dt due to R_3')
axis([0,10,-100000,100000])

saveas(energy_deriv,'energy_deriv','png')