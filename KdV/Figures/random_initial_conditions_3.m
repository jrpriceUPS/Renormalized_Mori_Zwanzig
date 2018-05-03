%produces coefficients, finds scaling law behavior for them, and plots the
%results to demonstrate how well they fit!

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

%use N values that will yield valid ROMs for all chosen epsilon values
N_list = 40:2:60;
epsilon = fliplr(0.065:0.005:0.1);
num_runs = 10;

t2 = zeros(length(N_list),length(epsilon),num_runs);
t4 = zeros(length(N_list),length(epsilon),num_runs);

t2_form = zeros(3,num_runs);
t4_form = zeros(3,num_runs);

a_list = zeros(num_runs,1);
b_list = zeros(num_runs,1);

for k = 1:num_runs
    
    a = rand;
    a_list(k) = a;
    b = rand*(1-a);
    b_list(k) = b;
    
    %run exact solution to time 10 and use to find t2 and t4 coefficients for
    %each ROM
    for j = 1:length(epsilon)
        epsilon(j)
        simulation_params.epsilon = epsilon(j);  %coefficient on linear term in KdV
        simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
        simulation_params.dt = 1e-3;      %timestep
        simulation_params.endtime = 10;   %end of simulation
        simulation_params.howoften = 1;   %how often to save state vector
        simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
        simulation_params.tol = inf;    %tolerance for identifying instabilities
        simulation_params.N = 256;          %number of positive modes to simulate
        simulation_params.order = 4;
        simulation_params.initial_condition = @(x) a*sin(x) + b*sin(2*x) + (1-a-b)*sin(3*x);
        
        %full model with no approximations
        simulation_params.initialization = @(x) full_init_KdV(x);
        
        [t_list,u_list] = PDE_solve(simulation_params);
        
        simulation_params = full_init(simulation_params);
        
        [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);
        coeffs_list = no_k_dependence_coeffs(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N_list,0,10,0);
        
        t2(:,j,k) = coeffs_list(:,2);
        t4(:,j,k) = coeffs_list(:,4);
        [coeffs_list(:,2) coeffs_list(:,4)]
        
    end
    
    %compute the scaling law using log-log least squares method
    t2_form(:,k) = find_scaling_law(squeeze(t2(:,:,k)).',epsilon,N_list);
    t4_form(:,k) = find_scaling_law(squeeze(t4(:,:,k)).',epsilon,N_list);
    
end

for k = 1:num_runs
    
    a = rand;
    a_list(k) = a;
    
    %run exact solution to time 10 and use to find t2 and t4 coefficients for
    %each ROM
    for j = 1:length(epsilon)
        epsilon(j)
        simulation_params.epsilon = epsilon(j);  %coefficient on linear term in KdV
        simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
        simulation_params.dt = 1e-3;      %timestep
        simulation_params.endtime = 10;   %end of simulation
        simulation_params.howoften = 1;   %how often to save state vector
        simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
        simulation_params.tol = inf;    %tolerance for identifying instabilities
        simulation_params.N = 256;          %number of positive modes to simulate
        simulation_params.order = 4;
        simulation_params.initial_condition = @(x) a*sin(x) + (1-a)*sin(2*x);
        
        %full model with no approximations
        simulation_params.initialization = @(x) full_init_KdV(x);
        
        [t_list,u_list] = PDE_solve(simulation_params);
        
        simulation_params = full_init_KdV(simulation_params);
        
        [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);
        coeffs_list = no_k_dependence_coeffs(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N_list,0,10,0);
        
        t2(:,j,k) = coeffs_list(:,2);
        t4(:,j,k) = coeffs_list(:,4);
        
    end
    
    %compute the scaling law using log-log least squares method
    t2_form2(:,k) = find_scaling_law(squeeze(t2(:,:,k)).',epsilon,N_list);
    t4_form2(:,k) = find_scaling_law(squeeze(t4(:,:,k)).',epsilon,N_list);
    
end





figure(1)
subplot(2,1,1)
plot(-t2_form(1,:),'.','markersize',20)
hold on
plot(-t2_form2(1,:),'r.','markersize',20)
legend('three modes activated','two modes activated')
xlabel('run','fontsize',16)
ylabel('t^2 prefactor','fontsize',16)

subplot(2,1,2)
plot(-t4_form(1,:),'.','markersize',20)
hold on
plot(-t4_form2(1,:),'r.','markersize',20)
legend('three modes activated','two modes activated')
xlabel('run','fontsize',16)
ylabel('t^4 prefactor','fontsize',16)
saveas(gcf,'prefactor3','png')


figure(2)
subplot(2,1,1)
plot(t2_form(2,:),'.','markersize',20)
hold on
plot(t2_form2(2,:),'r.','markersize',20)
legend('three modes activated','two modes activated')
xlabel('run','fontsize',16)
ylabel('t^2 N exponent','fontsize',16)

subplot(2,1,2)
plot(t4_form(2,:),'.','markersize',20)
hold on
plot(t4_form2(2,:),'r.','markersize',20)
legend('three modes activated','two modes activated')
xlabel('run','fontsize',16)
ylabel('t^4 N exponent','fontsize',16)
saveas(gcf,'N_exp3','png')


figure(3)
subplot(2,1,1)
plot(t2_form(3,:),'.','markersize',20)
hold on
plot(t2_form2(3,:),'r.','markersize',20)
legend('three modes activated','two modes activated')
xlabel('run','fontsize',16)
ylabel('t^2 epsilon exponent','fontsize',16)

subplot(2,1,2)
plot(t4_form(3,:),'.','markersize',20)
hold on
plot(t4_form2(3,:),'r.','markersize',20)
legend('three modes activated','two modes activated')
xlabel('run','fontsize',16)
ylabel('t^4 epsilon exponent','fontsize',16)
saveas(gcf,'eps_exp3','png')

