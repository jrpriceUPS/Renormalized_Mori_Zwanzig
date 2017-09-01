%produces coefficients, finds scaling law behavior for them, and plots the
%results to demonstrate how well they fit!

%clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

%use N values that will yield valid ROMs for all chosen epsilon values
N_list = 8:2:26;
epsilon = 0.1;

%run exact solution to time 10 and use to find t2 and t4 coefficients for
%each ROM


simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
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
    
    clear('simulation_params')
    
    simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
    simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
    simulation_params.dt = 1e-3;      %timestep
    simulation_params.endtime = 10;   %end of simulation
    simulation_params.howoften = 1;   %how often to save state vector
    simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
    simulation_params.tol = inf;    %tolerance for identifying instabilities
    simulation_params.N = N;          %number of positive modes to simulate
    simulation_params.initial_condition = @(x) sin(x);
    simulation_params.name = 'complete';
    
    simulation_params.M = 3*N;
    simulation_params.F_modes = [1:N,2*N:4*N+2,5*N+2:6*N];
    simulation_params.G_modes = N+1:5*N+1;
    simulation_params.k = [0:3*N-1,-3*N:-1].';
    [t_standard,u_standard] = KdV_solve(simulation_params);
    
    clear('simulation_params')
    
    simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
    simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
    simulation_params.dt = 1e-3;      %timestep
    simulation_params.endtime = 10;   %end of simulation
    simulation_params.howoften = 1;   %how often to save state vector
    simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
    simulation_params.tol = inf;    %tolerance for identifying instabilities
    simulation_params.N = N;          %number of positive modes to simulate
    simulation_params.initial_condition = @(x) sin(x);
    simulation_params.name = 'complete';
    
    simulation_params.M = 128;
    simulation_params.F_modes = [1:N,256-4*N:256-2*N+2,256-N+2:256];
    simulation_params.G_modes = N+1:256-N+1;
    simulation_params.k = [0:128-1,-128:-1].';
    [t_larger,u_larger] = KdV_solve(simulation_params);
    
    
    
    [x,u_real] = make_real_space(u_list(1:N,:),N);
    [~,u_standard_real] = make_real_space(u_standard,N);
    [~,u_larger_real] = make_real_space(u_larger,N);
    
    err_standard = sum((u_real(:,1:length(t_standard))-u_standard_real).^2,1)./sum(u_real(:,1:length(t_standard)).^2,1);
    err_larger = sum((u_larger_real-u_real(:,1:length(t_larger))).^2,1)./sum(u_real(:,1:length(t_larger)).^2,1);
    
    error = figure(2);
    set(gca,'FontSize',16)
    
    hold off
    plot(t_standard,err_standard,'r')
    hold on
    plot(t_larger,err_larger,'k')
    title(sprintf('Relative error of size N = %i models',N))
    xlabel('time')
    ylabel('relative global error')
    legend('M = 3N','M = 128','location','northwest')
    saveas(error,sprintf('M_error%i',N),'png')
    
end