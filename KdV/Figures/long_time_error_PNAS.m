%a script to produce plots of the error in a long time integration for a
%fixed epsilon and compare them against the markov model

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

N_list = 8:4:32;
epsilon = 0.1;
endtime = 100;
howoften = 100;

%find the exact solution
simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = endtime;   %end of simulation
simulation_params.howoften = howoften;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 256;          %number of positive modes to simulate
simulation_params.name = 'full';  %full simulation

simulation_params.initial_condition = @(x) sin(x);

[t_list,u_list] = KdV_solve(simulation_params);

u_list_markov = cell(length(N_list),1);
u_list4 = cell(length(N_list),1);

err_list_markov = zeros(length(N_list),1);
err_list4 = zeros(length(N_list),1);

for i = 1:length(N_list)
    
    N_list(i)
    
    simulation_params.name = 'full';
    simulation_params.N = N_list(i);
    [t_markov,u_markov] = KdV_solve(simulation_params);
    
    simulation_params.name = 'complete';
    simulation_params.order = 4;
    simulation_params.N = N_list(i);
    [t4,u4] = KdV_solve(simulation_params);
    
    
    u_list_markov{i} = u_markov;
    u_list4{i} = u4;
    
    
    figure(2)
    [x,u_real] = make_real_space(u_list(1:N_list(i),:),N_list(i));
    [~,u_markov_real] = make_real_space(u_markov,N_list(i));
    [~,u_4_real] = make_real_space(u4,N_list(i));
    
    err_markov = sum((u_real(:,1:length(t_markov))-u_markov_real).^2,1)./sum(u_real(:,1:length(t_markov)).^2,1);
    err4 = sum((u_4_real-u_real(:,1:length(t4))).^2,1)./sum(u_real(:,1:length(t4)).^2,1);

    
    if length(t_markov) == length(t_list)
        err_list_markov(i) = err_markov(end);
    else
        err_list_markov(i) = NaN;
    end
    
    if length(t4) == length(t_list)
        err_list4(i) = err4(end);
    else
        err_list4(i) = NaN;
    end
        
    
end

figure
set(gca,'FontSize',16)
semilogy(N_list,err_list_markov,'r.','markersize',20)
hold on
semilogy(N_list,err_list4,'k.','markersize',20)
legend('Markov','4th Order ROM','location','southwest')
xlabel('Resolution')
ylabel('Relative error at t = 100 (log plot)')
saveas(gcf,'logerr','png')





close all