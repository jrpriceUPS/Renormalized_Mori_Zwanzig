clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

sim_endtime = 1000;
N_full = 48;
N_list = 4:2:24;
N = 14;
tau_list = [0.0;0.2;0.4;0.6;0.8;1.0];

% Load calculated coefficients
% This data is generated using c_vs_N_test_coeff_data.m
load(['c4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])


for i = 1:length(tau_list)
    
    tau = tau_list(i);
    
    % Find coefficients corresponding to tau and N
    ind = find(tau == tau_list); 
    c4_op_tau = c4_op{1,ind};
    ind_N = find(N == N_list); 
    c4_op_N_tau = c4_op_tau(:,ind_N);
        
    % Load t and u data for specified N and tau
    % This data is generated using c_vs_N_tests_sims.m
    load(['t_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
    load(['u_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
    
    % Check that the sizes match
    if size(u,1) ~= N
        break
    end
    
    % Calculate the energy in the N modes of the ROM - use as a
    % verification step
    energy = get_3D_energy(u,N);
    
    % Calculate the change in energy of each of the memory terms
    [R0_test,R1_test,R2_test,R3_test,R4_test,RTS_test,RT_test] = ROM_dE_memory(c4_op_N_tau,N,t,u);
    
    figure()
    
    hold on
    c_list = (4/5)*hsv(5);
    plot(log(t),R0_test,'DisplayName','\DeltaE^0','Color',c_list(1,:))
    plot(log(t),R1_test,'DisplayName','\DeltaE^1','Color',c_list(2,:))
    plot(log(t),R2_test,'DisplayName','\DeltaE^2','Color',c_list(3,:))
    plot(log(t),R3_test,'DisplayName','\DeltaE^3','Color',c_list(4,:))
    plot(log(t),R4_test,'DisplayName','\DeltaE^4','Color',c_list(5,:))
    box on
    xlim([log(0.1),log(sim_endtime)])
    xlabel('log(t)')
    ylabel('\DeltaE^n')
    hold off
    legend show
    
    saveas(gcf,sprintf('Euler_memory_energy_%i_M_%i_N_%i_tau_%i_c4_test',sim_endtime,N_full,N,tau),'png')
    
    figure()
    wstart = 2.5;
    xmin = t(find(t >= wstart));
    wend = 20;
    xmax = t(find(t >= wend));
    hold on
    c_list = (4/5)*hsv(5);
    plot(log(t(find(t >= wstart & t <= wend))),R0_test(find(t >= wstart & t <= wend)),'DisplayName','\DeltaE^0','Color',c_list(1,:))
    plot(log(t(find(t >= wstart & t <= wend))),R3_test(find(t >= wstart & t <= wend)),'DisplayName','\DeltaE^3','Color',c_list(4,:))
    plot(log(t(find(t >= wstart & t <= wend))),R4_test(find(t >= wstart & t <= wend)),'DisplayName','\DeltaE^4','Color',c_list(5,:))
    box on
    xlim([log(xmin(1)),log(xmax(1))])
    xlabel('log(t)')
    ylabel('\DeltaE^n')
    hold off
    legend('location','southwest')
    legend show
    
    saveas(gcf,sprintf('Euler_memory_energy_%i_M_%i_N_%i_tau_%i_c4_short_times_test',sim_endtime,N_full,N,tau),'png')
    
    
    % Save the results into the directory
    save(['R0_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'R0_test');
    save(['R1_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'R1_test');
    save(['R2_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'R2_test');
    save(['R3_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'R3_test');
    save(['R4_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'R4_test');
    save(['RT_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'RT_test');

end
