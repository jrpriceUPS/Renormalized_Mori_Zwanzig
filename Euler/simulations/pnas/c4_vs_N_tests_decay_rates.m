clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

N_full = 48;
sim_endtime = 1000;
alpha = 1;
N_list = [6:2:14];
tau_list = 1.0;

initial_decay_time{1,length(N_list)} = [];
initial_decay_rate{1,length(N_list)} = [];
second_decay_rate{1,length(N_list)} = [];

for i = 1:length(N_list)
    
    N = N_list(i);
    initial_decay_time{i} = zeros(length(tau_list),1);
    initial_decay_rate{i} = zeros(length(tau_list),1);
    second_decay_rate{i} = zeros(length(tau_list),1);
    
    for j = 1:length(tau_list)
        
        tau = tau_list(j);
        
        % This data is generated using c_vs_N_tests_sims.m
        load(['t_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
        load(['u_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
        
        % Check that the sizes match
        if size(u,1) ~= N
            break
        end
        
        % Calculate the energy in the N modes of the ROM
        energy = get_3D_energy(u,N);
        norm_energy = 1 - energy./energy(1);
        
        % Find 10%, 50%, 90% and 99.5% energy locations
        ind_10 = max(find(norm_energy <= 0.10));
        ind_50 = min(find(norm_energy >= 0.50));
        ind_90 = max(find(norm_energy <= 0.90));
        ind_995 = min(find(norm_energy >= 0.995));
        
        initial_decay_time{i}(j,1) = t(ind_10);
        
        int_decay_fit = polyfit(log(t(ind_50:ind_90)),log(energy(ind_50:ind_90)), 1);
        initial_decay_rate{i}(j,1) = int_decay_fit(1);
        
        sec_decay_fit = polyfit(log(t(ind_995:end)),log(energy(ind_995:end)), 1);
        second_decay_rate{i}(j,1) = sec_decay_fit(1);
        
    end
    
end

