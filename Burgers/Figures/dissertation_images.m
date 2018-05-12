% creates images for my dissertation

clear all;close all

% set of resolutions to present
N_list = 6:2:12;
endtime = 1000;

% compute the data!
leg = {'Exact','n = 1, constant','n = 2, constant','n = 3, constant','n = 4, constant','n = 1, decaying','n = 2, decaying','n = 3, decaying','n = 4, decaying'};
[times_array,energies_array,error_array] = generate_comparisons(N_list,endtime);


% Create energy plots for constant and algebraically decaying 
% renormalization coefficient ROMs
for i = 1:4
    
    N = N_list(i);
    
    t_list = times_array{i}.exact;
    tc1B = times_array{i}.c1B;
    tc2B = times_array{i}.c2B;
    tc3B = times_array{i}.c3B;
    tc4B = times_array{i}.c4B;
    
    tc1KdV = times_array{i}.c1KdV;
    tc2KdV = times_array{i}.c2KdV;
    tc3KdV = times_array{i}.c3KdV;
    tc4KdV = times_array{i}.c4KdV;
    
    energy_exact = energies_array{i}.exact;
    energyc1B = energies_array{i}.c1B;
    energyc2B = energies_array{i}.c2B;
    energyc3B = energies_array{i}.c3B;
    energyc4B = energies_array{i}.c4B;
    
    energyc1KdV = energies_array{i}.c1KdV;
    energyc2KdV = energies_array{i}.c2KdV;
    energyc3KdV = energies_array{i}.c3KdV;
    energyc4KdV = energies_array{i}.c4KdV;
    
    errc1B = error_array{i}.c1B;
    errc2B = error_array{i}.c2B;
    errc3B = error_array{i}.c3B;
    errc4B = error_array{i}.c4B;
    
    figure(1);
    subplot(2,2,i)
    hold off
    plot(log(t_list),log(energy_exact),'linewidth',2)
    hold on
    plot(log(tc1B),log(energyc1B),'r')
    plot(log(tc2B),log(energyc2B),'k')
    plot(log(tc3B),log(energyc3B),'c')
    plot(log(tc4B),log(energyc4B),'m')
    
    plot(log(tc1KdV),log(energyc1KdV),'r--','linewidth',1.2)
    plot(log(tc2KdV),log(energyc2KdV),'k--','linewidth',1.2)
    plot(log(tc3KdV),log(energyc3KdV),'c--','linewidth',1.2)
    plot(log(tc4KdV),log(energyc4KdV),'m--','linewidth',1.2)
    lgd = legend(leg{:},'location','southwest');
    lgd.FontSize = 6;
    axis([log(0.1),log(1000),-15,0])
    
    title(sprintf('N = %i',N),'fontsize',16)
    xlabel('log(t)')
    ylabel('log(energy)')
    
    
end

saveas(gcf,'energy_burgers','png')



% Create error for constant renormalization coefficient ROMs
for i = 1:4
    
     N = N_list(i);
    
    t_list = times_array{i}.exact;
    tc1B = times_array{i}.c1B;
    tc2B = times_array{i}.c2B;
    tc3B = times_array{i}.c3B;
    tc4B = times_array{i}.c4B;
    
    errc1B = error_array{i}.c1B;
    errc2B = error_array{i}.c2B;
    errc3B = error_array{i}.c3B;
    errc4B = error_array{i}.c4B;
    
    figure(2);
    subplot(2,2,i)
    hold off
    plot(tc1B,errc1B,'r','linewidth',1.5)
    hold on
    plot(tc2B,errc2B,'k','linewidth',1.5)
    plot(tc3B,errc3B,'c','linewidth',1.5)
    plot(tc4B,errc4B,'m','linewidth',1.5)
    axis([0,endtime,0,3])
    legend(leg{2:5},'location','northeast')
    
    title(sprintf('N = %i',N),'fontsize',16)
    xlabel('t')
    ylabel('error')
    
end

saveas(gcf,'error_burgers','png')