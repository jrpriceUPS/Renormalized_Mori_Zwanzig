%A script to plot how the k-dependent optimal coefficients (by least
%squares fit) compare with a shifting window for the fitting. This is
%compared to the formula for coefficients. This is useful for visually
%seeing how the "optimal coefficients" are largely static and captured well
%by our formulae

close all
width = 5;
spacing = 0.1;

coeffs_list = k_dependent_coefficients3(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,width,spacing,0);
coeffs_list2 = no_k_dependence_coeffs2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,0);

for k = 1:9;
    
    for i = 1:(10-width)/spacing+1
        subplot(2,1,1)
        hold off
        plot(1:N_list(k)-1,coeffs_list(k,1:N_list(k)-1,i))
        hold on
        plot([1,N_list(k)-1],ones(2,1)*coeffs_list2(k,1),'r')
        title(sprintf('N = %i',N_list(k)),'fontsize',16)
        
        subplot(2,1,2)
        hold off
        plot(1:N_list(k)-1,coeffs_list(k,1*(N_list(k)-1)+1:2*(N_list(k)-1),i))
        hold on
        plot([1,N_list(k)-1],ones(2,1)*coeffs_list2(k,2),'r')
        title(sprintf('N = %i',N_list(k)),'fontsize',16)
        
        
        
        pause(0.1)
    end
    pause
end