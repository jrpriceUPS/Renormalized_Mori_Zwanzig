%a script to produce plots of the energy change for N = 20

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

N = 20;
epsilon = 0.1;
endtime = 10;
howoften = 10;

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

simulation_params.name = 'full';
simulation_params.N = N;
[t_markov,u_markov] = KdV_solve(simulation_params);

simulation_params.name = 'complete';
simulation_params.order = 4;
simulation_params.N = N;
[t4,u4] = KdV_solve(simulation_params);

%parameter details for non-renormalized ROM simulation
simulation_params.N = N;               %number of positive modes to simulate
simulation_params.name = 'complete';    %complete ROM
simulation_params.order = 4;            %use fourth order ROM
simulation_params.time_dependence = 1;  %include time dependence!
simulation_params.coeffs = ones(4,1);   %no renormalization
simulation_params.dt = 1e-5;


[t_blowup,u_blowup] = KdV_solve(simulation_params);
clear('simulation_params')


figure(1)
set(gca,'FontSize',16)
hold off
plot(t_list,get_energy(u_list,N),'linewidth',2)
hold on
plot(t_markov,get_energy(u_markov,N),'r','linewidth',2)
plot(t4,get_energy(u4,N),'k','linewidth',2)
plot(t_blowup,get_energy(u_blowup,N),'c','linewidth',2)
ax = axis;
axis([0,10,0.4982,0.5002]);
legend('Exact','Markov','4th order ROM','Non-renormalized 4th order ROM','location','southwest')
xlabel('time','fontsize',16)
ylabel('Mass in resolved modes','fontsize',16)
title(sprintf('N = %i',N),'fontsize',16)
saveas(gcf,sprintf('energy_evo%i',N),'png')

figure(2)
set(gca,'FontSize',16)
hold off
plot(t_list,get_energy(u_list,N),'linewidth',2)
hold on
plot(t_markov,get_energy(u_markov,N),'r','linewidth',2)
plot(t_blowup,get_energy(u_blowup,N),'c','linewidth',2)
ax = axis;
axis([0,10,0.4994,0.5002]);
legend('Exact','Markov','Non-renormalized complete ROM','location','southwest')
xlabel('time','fontsize',16)
ylabel('Mass in resolved modes','fontsize',16)
title(sprintf('N = %i',N),'fontsize',16)
saveas(gcf,'blowup','png')

figure(3)
set(gca,'FontSize',16)
hold off
plot(t_list,get_energy(u_list,8),'linewidth',2)
axis([0,10,0.34,0.54])
legend('Mass in first 8 Fourier modes')
saveas(gcf,'kdv_energy','png')

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = endtime;   %end of simulation
simulation_params.howoften = howoften;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 8;          %number of positive modes to simulate
simulation_params.name = 'full';  %full simulation
[t_markov8,u_markov8] = KdV_solve(simulation_params);

figure(4)
set(gca,'FontSize',16)
hold off
plot(t_list,get_energy(u_list,8),'linewidth',2)
hold on
plot(t_markov8,get_energy(u_markov8,8),'r','linewidth',2)
axis([0,10,0.34,0.54])
legend('Mass in first 8 Fourier modes (exact)','Mass in first 8 Fourier modes (Markov)')
saveas(gcf,'markov_energy','png')




close all