%generates plot showing that the non-renormalized models blow up!

clear all;close all;
N = 20;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

%parameter details for non-renormalized ROM simulation
simulation_params.epsilon = 0.1;       %coefficient on linear term in KdV
simulation_params.alpha = 1;            %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;            %timestep
simulation_params.endtime = 10;         %end of simulation
simulation_params.howoften = 1;         %how often to save state vector
simulation_params.blowup = 1;           %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;            %tolerance for identifying instabilities
simulation_params.N = N;               %number of positive modes to simulate
simulation_params.initialization = @(x) complete_init_KdV(x);    %complete ROM
simulation_params.order = 4;            %use fourth order ROM
simulation_params.time_dependence = 1;  %include time dependence!
simulation_params.coeffs = [1;-1/2;1/6;-1/24];   %no renormalization

[t_list,u_list] = PDE_solve(simulation_params);
clear('simulation_params')


%parameter details for non-renormalized ROM simulation
simulation_params.epsilon = 0.1;       %coefficient on linear term in KdV
simulation_params.alpha = 1;            %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;            %timestep
simulation_params.endtime = 10;         %end of simulation
simulation_params.howoften = 1;         %how often to save state vector
simulation_params.blowup = 1;           %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;            %tolerance for identifying instabilities
simulation_params.N = N;               %number of positive modes to simulate
simulation_params.initialization = @(x) complete_init_KdV(x);    %complete ROM
simulation_params.order = 2;            %use fourth order ROM
simulation_params.time_dependence = 1;  %include time dependence!
simulation_params.coeffs = [1;-1/2];   %no renormalization

[t_list2,u_list2] = PDE_solve(simulation_params);
clear('simulation_params')



%parameter details for exact simulation
simulation_params.epsilon = 0.1;       %coefficient on linear term in KdV
simulation_params.alpha = 1;            %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;            %timestep
simulation_params.endtime = 10;         %end of simulation
simulation_params.howoften = 1;         %how often to save state vector
simulation_params.blowup = 1;           %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;            %tolerance for identifying instabilities
simulation_params.N = 128;              %number of positive modes to simulate
simulation_params.initialization = @(x) full_init_KdV(x);        %complete ROM

[t_exact,u_exact] = PDE_solve(simulation_params);
clear('simulation_params')

plot(t_exact,get_energy(u_exact,N))
ax = axis;
hold on
plot(t_list2,get_energy(u_list2,N),'r')
plot(t_list,get_energy(u_list,N),'k')
axis([0,10,ax(3),ax(4)])
title(sprintf('N = %i',N),'fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy in resolved modes','fontsize',16)
legend('exact','non-renormalized 2nd order ROM','non-renormalized 4th order ROM','location','southeast')
saveas(gcf,'blowup','png')