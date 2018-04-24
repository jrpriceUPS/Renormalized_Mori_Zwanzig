clear all;close all;

addpath ../renormalization_test_Apr17
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N_list = 10:2:24;
alpha = 1;
endtime = 100;
if ~(exist(sprintf('u_list%i.mat',endtime),'file') == 2)
    
    [t_list,u_list] = upwind_burgers(1,10000,endtime,1e-4,100);
    save(sprintf('u_list%i',endtime),'u_list');
    save(sprintf('t_list%i',endtime),'t_list');
     
else
     load(sprintf('u_list%i',endtime));
     load(sprintf('t_list%i',endtime));
end

leg = {'exact','c1B','c3B13','c2B','c3B','c4B','c1KdV','c3KdV13','c2KdV','c3KdV','c4KdV'};


for i = 1:length(N_list)
    
    N = N_list(i);
    
    simulation_params.N = N;
    simulation_params.alpha = alpha;
    simulation_params.endtime = endtime;
    simulation_params.print_time = 1;
    
    % Burgers-style 1 ROM
    type = 'c1B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc1B,uc1B] = PDE_solve(simulation_params);
    
    % Burgers-style 1+3 ROM
    type = 'c3B13';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc3B13,uc3B13] = PDE_solve(simulation_params);
    
    % Burgers-style 1+2 ROM
    type = 'c2B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc2B,uc2B] = PDE_solve(simulation_params);
    
    % Burgers-style 1+2+3 ROM
    type = 'c3B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc3B,uc3B] = PDE_solve(simulation_params);
    
    % Burgers-style 1+2+3+4 ROM
    type = 'c4B';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 1;
    [tc4B,uc4B] = PDE_solve(simulation_params);
    
    
    
    % KdV-style 1 ROM
    type = 'c1KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc1KdV,uc1KdV] = PDE_solve(simulation_params);
    
    % KdV-style 1+3 ROM
    type = 'c3KdV13';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc3KdV13,uc3KdV13] = PDE_solve(simulation_params);
    
    % KdV-style 1+2 ROM
    type = 'c2KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc2KdV,uc2KdV] = PDE_solve(simulation_params);
    
    % KdV-style 1+2+3 ROM
    type = 'c3KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc3KdV,uc3KdV] = PDE_solve(simulation_params);
    
    % KdV-style 1+2+3+4 ROM
    type = 'c4KdV';
    simulation_params.time_dependence = 1;
    simulation_params.degree = 4;
    simulation_params.coeffs = scaling_law(N,type);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    simulation_params.time_dependence = 0;
    [tc4KdV,uc4KdV] = PDE_solve(simulation_params);
    
    
    
    
    
    
    energy_exact = get_energy(u_list,N);
    energyc1B = get_energy(uc1B,N);
    energyc3B13 = get_energy(uc3B13,N);
    energyc2B = get_energy(uc2B,N);
    energyc3B = get_energy(uc3B,N);
    energyc4B = get_energy(uc4B,N);
    energyc1KdV = get_energy(uc1KdV,N);
    energyc3KdV13 = get_energy(uc3KdV13,N);
    energyc2KdV = get_energy(uc2KdV,N);
    energyc3KdV = get_energy(uc3KdV,N);
    energyc4KdV = get_energy(uc4KdV,N);
    
    figure(1)
    hold off
    plot(log(t_list),log(energy_exact),'linewidth',2)
    hold on
    plot(log(tc1B),log(energyc1B),'r')
    plot(log(tc3B13),log(energyc3B13),'g')
    plot(log(tc2B),log(energyc2B),'k')
    plot(log(tc3B),log(energyc3B),'c')
    plot(log(tc4B),log(energyc4B),'m')
    
    plot(log(tc1KdV),log(energyc1KdV),'r--','linewidth',1.2)
    plot(log(tc3KdV13),log(energyc3KdV13),'g--','linewidth',1.2)
    plot(log(tc2KdV),log(energyc2KdV),'k--','linewidth',1.2)
    plot(log(tc3KdV),log(energyc3KdV),'c--','linewidth',1.2)
    plot(log(tc4KdV),log(energyc4KdV),'m--','linewidth',1.2)
    legend(leg{:},'location','southwest')
    
    title(sprintf('N = %i',N),'fontsize',16)
    xlabel('log(t)')
    ylabel('log(energy)')
    saveas(gcf,sprintf('energy%i',N),'png')
    close
    
end