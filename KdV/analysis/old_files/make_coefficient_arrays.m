%produces arrays of k-dependent coefficients for easier analysis (probably
%won't be used much anymore)
clear all;close all;

%load locally save data
load u_list
load t_list
load simulation_params
load N_list

%compute energy derivatives and find k-dependent coefficients
[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);
coeffs_list = k_dependent_coefficients2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N_list,0,10,0);

alpha_array = zeros(length(N_list),N_list(end)-1);
beta_array = zeros(length(N_list),N_list(end)-1);
gamma_array = zeros(length(N_list),N_list(end)-1);
eta_array = zeros(length(N_list),N_list(end)-1);

for i = 1:length(N_list)
    
    N = N_list(i);
    
    alpha = coeffs_list(i,1:(N-1)); %t-model
    beta = coeffs_list(i,(N-1)+1:2*(N-1)); %t2-model
    gamma = coeffs_list(i,2*(N-1)+1:3*(N-1)); %t3-model
    eta = coeffs_list(i,3*(N-1)+1:4*(N-1)); %t4-model
    
    alpha_array(i,1:N-1) = alpha;
    beta_array(i,1:N-1) = beta;
    gamma_array(i,1:N-1) = gamma;
    eta_array(i,1:N-1) = eta;
    
end