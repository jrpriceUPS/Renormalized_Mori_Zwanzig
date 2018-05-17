%a script to produce plots of the "locking in" of coefficients

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

%use N values that will yield valid ROMs for all chosen epsilon values
N_list = 32:2:60;
epsilon = fliplr(0.065:0.005:0.1);

endtime = 10;
skip = 0.1;
t = 2:skip:endtime;

c1 = zeros(length(N_list),length(epsilon),length(t));
c2 = zeros(length(N_list),length(epsilon),length(t));
c3 = zeros(length(N_list),length(epsilon),length(t));
c4 = zeros(length(N_list),length(epsilon),length(t));

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
    simulation_params.initial_condition = @(x) sin(x);
    
    %full model with no approximations
    simulation_params.initialization = @(x) full_init_KdV(x);
    
    [t_list,u_list] = PDE_solve(simulation_params);
    
    simulation_params = full_init_KdV(simulation_params);
    
    [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);
    coeffs_list = no_k_dependence_coeffs(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N_list,0,.1,0);
    
    c1(:,j,:) = coeffs_list(:,1,:);
    c2(:,j,:) = coeffs_list(:,2,:);
    c3(:,j,:) = coeffs_list(:,3,:);
    c4(:,j,:) = coeffs_list(:,4,:);
    
end

t1_form = zeros(3,length(t));
t2_form = zeros(3,length(t));
t3_form = zeros(3,length(t));
t4_form = zeros(3,length(t));

for i = 1:length(t)
    t1_form(:,i) = find_scaling_law(squeeze(c1(:,:,i)),epsilon,N_list);
    t2_form(:,i) = find_scaling_law(squeeze(c2(:,:,i)),epsilon,N_list);
    t3_form(:,i) = find_scaling_law(squeeze(c3(:,:,i)),epsilon,N_list);
    t4_form(:,i) = find_scaling_law(squeeze(c4(:,:,i)),epsilon,N_list);
end

figure(1)
subplot(3,1,1)
plot(t,t1_form(1,:))
title('t-model coefficient')
xlabel('time')
ylabel('Prefactor')
subplot(3,1,2)
plot(t,t1_form(2,:))
xlabel('time')
ylabel('Exponent on N')
subplot(3,1,3)
plot(t,t1_form(3,:))
xlabel('Edge of window [2,t]')
ylabel('Exponent on epsilon')
saveas(gcf,'t1_coeffs','png')


figure(2)
subplot(3,1,1)
plot(t,t2_form(1,:))
title('t^2-model coefficient')
xlabel('time')
ylabel('Prefactor')
subplot(3,1,2)
plot(t,t2_form(2,:))
xlabel('time')
ylabel('Exponent on N')
subplot(3,1,3)
plot(t,t2_form(3,:))
xlabel('Edge of window [2,t]')
ylabel('Exponent on epsilon')
saveas(gcf,'t2_coeffs','png')


figure(3)
subplot(3,1,1)
plot(t,t3_form(1,:))
title('t^3-model coefficient')
xlabel('time')
ylabel('Prefactor')
subplot(3,1,2)
plot(t,t3_form(2,:))
xlabel('time')
ylabel('Exponent on N')
subplot(3,1,3)
plot(t,t3_form(3,:))
xlabel('Edge of window [2,t]')
ylabel('Exponent on epsilon')
saveas(gcf,'t3_coeffs','png')


figure(4)
subplot(3,1,1)
plot(t,t4_form(1,:))
title('t^4-model coefficient')
xlabel('time')
ylabel('Prefactor')
subplot(3,1,2)
plot(t,t4_form(2,:))
xlabel('time')
ylabel('Exponent on N')
subplot(3,1,3)
plot(t,t4_form(3,:))
xlabel('Edge of window [2,t]')
ylabel('Exponent on epsilon')
saveas(gcf,'t4_coeffs','png')