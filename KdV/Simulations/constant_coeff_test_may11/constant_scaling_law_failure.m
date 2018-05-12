%clear all;close all;


addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N = 20;
epsilon = 0.1;

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 256;          %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);

%full model with no approximations
simulation_params.initialization = @(x) full_init_KdV(x);

[t_list,u_list] = PDE_solve(simulation_params);

simulation_params = full_init_KdV(simulation_params);

[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N);
coeffs_list = constant_coeffs(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N,0,10,1);
coeffs_list2 = coeffs_list;
coeffs_list2(1) = 0;
coeffs_list2(3) = 0;

clear('simulation_params')




simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 10;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 256;          %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);

%full model with no approximations
simulation_params.initialization = @(x) full_init_KdV(x);

[t_list,u_list] = PDE_solve(simulation_params);


simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 10;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.initial_condition = @(x) sin(x);
simulation_params.initialization = @(x) complete_init_constant_KdV(x);
simulation_params.coeffs = coeffs_list;
simulation_params.order = 4;
simulation_params.N = N;
simulation_params.time_dependence = 1;
[t_all,u_all] = PDE_solve(simulation_params);

simulation_params.initialization = @(x) complete_init_constant_KdV(x);
simulation_params.coeffs = coeffs_list2;
simulation_params.order = 4;
simulation_params.N = N;
simulation_params.time_dependence = 1;
[t_24,u_24] = PDE_solve(simulation_params);


exact = get_energy(u_list,N);
all4 = get_energy(u_all,N);
just24 = get_energy(u_24,N);

figure(1)
hold off
plot(t_list,exact,'b')
hold on
plot(t_all,all4,'r')
plot(t_24,just24,'k')
title('N = 20','fontsize',16)
legend('exact','t-model through t^4-model','only t^2-model and t^4-model','location','northeast')
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
axis([0,100,0.499,0.51])
saveas(gcf,'constant_coeff_energy','png')



[x,u_real] = make_real_space(u_list(1:N,:),N);
[~,u_all_real] = make_real_space(u_all,N);
[~,u_24_real] = make_real_space(u_24,N);

err_all = sum((u_real(:,1:length(t_all))-u_all_real).^2,1)./sum(u_real(:,1:length(t_all)).^2,1);
err_24 = sum((u_24_real-u_real(:,1:length(t_24))).^2,1)./sum(u_real(:,1:length(t_24)).^2,1);

figure(2)
hold off
plot(t_all,err_all,'r')
plot(t_24,err_24,'k')
legend('t-model through t^4-model','only t^2-model and t^4-model','location','northeast')
xlabel('time','fontsize',16)
ylabel('relative error','fontsize',16)
title('N = 20','fontsize',16)

