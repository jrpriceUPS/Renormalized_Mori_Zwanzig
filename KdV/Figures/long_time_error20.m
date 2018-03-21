%a script to produce plots of the error in a long time integration for a
%fixed epsilon and compare them against the markov model

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

N = 20;
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

u_list_markov = cell(length(N),1);
u_list4 = cell(length(N),1);
u_list2 = cell(length(N),1);



simulation_params.name = 'full';
simulation_params.N = N;
[t_markov,u_markov] = KdV_solve(simulation_params);

simulation_params.name = 'complete';
simulation_params.order = 4;
simulation_params.N = N;
[t4,u4] = KdV_solve(simulation_params);

simulation_params.name = 'complete';
simulation_params.order = 2;
simulation_params.N = N;
[t2,u2] = KdV_solve(simulation_params);


figure(2)
[x,u_real] = make_real_space(u_list(1:N,:),N);
[~,u_markov_real] = make_real_space(u_markov,N);
[~,u_2_real] = make_real_space(u2,N);
[~,u_4_real] = make_real_space(u4,N);

err_markov = sum((u_real(:,1:length(t_markov))-u_markov_real).^2,1)./sum(u_real(:,1:length(t_markov)).^2,1);
err2 = sum((u_2_real-u_real(:,1:length(t2))).^2,1)./sum(u_real(:,1:length(t2)).^2,1);
err4 = sum((u_4_real-u_real(:,1:length(t4))).^2,1)./sum(u_real(:,1:length(t4)).^2,1);

hold off
plot(t_markov,err_markov,'k.')
hold on
plot(t2,err2,'k*')
plot(t4,err4,'ko')
legend('Markov','2nd order ROM','4th order ROM','location','northwest')
xlabel('time','fontsize',16)
ylabel('Real space relative global error','fontsize',16)
title(sprintf('N = %i',N),'fontsize',16)
saveas(gcf,sprintf('real_err%i',N),'png')



close all