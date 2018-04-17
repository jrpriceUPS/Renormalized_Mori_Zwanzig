%produces arrays of coefficients for easier analysis
clear all;close all;

load u_list
load t_list
load simulation_params
load N_list

[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow] = generate_deriv_data_4func2(t_list,u_list,simulation_params,N_list);
coeffs_list = k_dependent_coefficients3(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,0);


coeffs_list2 = no_k_dependence_coeffs2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,0);

alpha_array = zeros(length(N_list),N_list(end)-1);
beta_array = zeros(length(N_list),N_list(end)-1);

for i = 1:length(N_list)
    
    N = N_list(i);
    
    alpha = coeffs_list(i,1:(N-1)); %t-model
    beta = coeffs_list(i,(N-1)+1:2*(N-1)); %t2-model
    
    alpha_array(i,1:N-1) = alpha;
    beta_array(i,1:N-1) = beta;
    
end