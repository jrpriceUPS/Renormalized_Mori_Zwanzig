function coeffs_array = changing_optimal_k_BCH(t_list,u_list,simulation_params,width,spacing,N_list)

%A script to plot how the k-dependent optimal coefficients (by least
%squares fit) compare with a shifting window for the fitting. This is
%compared to the formula for coefficients. This is useful for visually
%seeing how the "optimal coefficients" are largely static and captured well
%by our formulae

close all

[u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow] = generate_deriv_data_BCH(t_list,u_list,simulation_params,N_list);

coeffs_list = k_dependent_coefficients3(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,width,spacing,0);
coeffs_array = no_k_dependence_coeffs2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,0);

pause
for k = 1:9;
    
    coeffs = coeffs_array(k,:);
    
    for i = 1:(max(t_list)-width)/spacing+1
        subplot(2,1,1)
        hold off
        plot(1:N_list(k)-1,coeffs_list(k,1:N_list(k)-1,i))
        hold on
        plot(1:N_list(k)-1,coeffs(1)*ones(1,N_list(k)-1),'r')
        title(sprintf('N = %i',N_list(k)),'fontsize',16)
        
        subplot(2,1,2)
        hold off
        plot(1:N_list(k)-1,coeffs_list(k,1*(N_list(k)-1)+1:2*(N_list(k)-1),i))
        hold on
        plot(1:N_list(k)-1,coeffs(2)*ones(1,N_list(k)-1),'r')
        title(sprintf('N = %i',N_list(k)),'fontsize',16)
        
        pause(0.1)
    end
    pause
end