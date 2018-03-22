%a function to compute a stable ROM for epsilon = 0.01 by extrapolating the
%scaling laws for the renormalization coefficients and compare to the exact
%answer

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

clear all;close all;
M = 512;
N = 256;
epsilon = 0.01;

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-4;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 100;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = M;          %number of positive modes to simulate
simulation_params.initial_condition = @(x) sin(x);
simulation_params.initialization = @(x) full_init_KdV(x);  %full simulation

[t_list,u_list] = PDE_solve(simulation_params);

if length(t_list) ~= length(0:simulation_params.dt*simulation_params.howoften:simulation_params.endtime)
    return
end


simulation_params.initialization = @(x) complete_init_KdV(x);
simulation_params.order = 4;
simulation_params.N = N;

[t_ROM,u_ROM] = PDE_solve(simulation_params);
save t_ROM t_ROM
save u_ROM u_ROM

simulation_params.initialization = @(x) full_init_KdV(x);
simulation_params.N = N;

[t_markov,u_markov] = PDE_solve(simulation_params);
save t_markov t_markov
save u_markov u_markov



energy = figure(1);
set(gca,'FontSize',16)
hold off
plot(t_list,get_energy(u_list,N));
hold on
plot(t_markov,get_energy(u_markov,N),'r');
plot(t_ROM,get_energy(u_ROM,N),'k');
title(sprintf('Mass in first N = %i modes',N))
xlabel('time')
ylabel('mass')
legend('Exact','Markov','Order 4 ROM','location','southwest')
saveas(energy,'extrap_energy','png')


[x,u_real] = make_real_space(u_list(1:N,:),N);
[~,u_markov_real] = make_real_space(u_markov,N);
[~,u_ROM_real] = make_real_space(u_ROM,N);

err_markov = sum((u_real(:,1:length(t_markov))-u_markov_real).^2,1)./sum(u_real(:,1:length(t_markov)).^2,1);
err_ROM = sum((u_ROM_real-u_real(:,1:length(t_ROM))).^2,1)./sum(u_real(:,1:length(t_ROM)).^2,1);

save u_list u_list
save t_list t_list

save u_markov u_markov
save t_markov t_markov

save u_ROM u_ROM
save t_ROM t_ROM

error = figure(2);
set(gca,'FontSize',16)

hold off
plot(t_markov,err_markov,'k:','linewidth',1.5)
hold on
plot(t_ROM,err_ROM,'k')
title(sprintf('Relative error of size N = %i models',N))
xlabel('time')
ylabel('relative global error')
legend('Markov','4th Order ROM','location','northwest')
saveas(error,'extrap_error','png')
close all