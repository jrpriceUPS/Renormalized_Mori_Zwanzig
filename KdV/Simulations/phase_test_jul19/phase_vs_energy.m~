%a function to compare the error in simulations using the coefficients from
%energy fits and mode fits

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

%clear all;close all;
N_list = 8:4:40;
epsilon = 0.1;

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 128;          %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);

%full model with no approximations
simulation_params.name = 'full';

[t_list,u_list] = KdV_solve(simulation_params);

for i = 1:length(N_list)
    
    N = N_list(i);
    
    simulation_params.name = 'complete';
    simulation_params.N = N;
    simulation_params.coeff_form = 'energy';
    
    [t_energy,u_energy] = KdV_solve(simulation_params);
    
    simulation_params.coeff_form = 'phase';
    
    [t_phase,u_phase] = KdV_solve(simulation_params);
    
    
    
    
    energy = figure(1);
    set(gca,'FontSize',16)
    hold off
    plot(t_list,get_energy(u_list,N));
    hold on
    plot(t_energy,get_energy(u_energy,N),'r');
    plot(t_phase,get_energy(u_phase,N),'k');
    title(sprintf('Mass in first N = %i modes',N))
    xlabel('time')
    ylabel('mass')
    legend('Exact','Energy fit ROM','Phase fit ROM','location','southwest')
    saveas(energy,sprintf('phase_energy%i',N),'png')
    
    
    [x,u_real] = make_real_space(u_list(1:N,:),N);
    [~,u_energy_real] = make_real_space(u_energy,N);
    [~,u_phase_real] = make_real_space(u_phase,N);
    
    err_energy = sum((u_real(:,1:length(t_energy))-u_energy_real).^2,1)./sum(u_real(:,1:length(t_energy)).^2,1);
    err_phase = sum((u_phase_real-u_real(:,1:length(t_phase))).^2,1)./sum(u_real(:,1:length(t_phase)).^2,1);
    
    
    error = figure(2);
    set(gca,'FontSize',16)
    
    hold off
    plot(t_energy,err_energy,'r')
    hold on
    plot(t_phase,err_phase,'k')
    title(sprintf('Relative error of size N = %i models',N))
    xlabel('time')
    ylabel('relative global error')
    legend('Energy fit ROM','Phase fit ROM','location','northwest')
    saveas(error,sprintf('phase_error%i',N),'png')
    
end