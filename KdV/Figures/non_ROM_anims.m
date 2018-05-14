% produces animations of the full solution and a comparison of the full
% solution to the Markov solution

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

N = 12;
epsilon = 0.1;

%exact solution

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 100;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 128;          %number of positive modes to simulate

%full model with no approximations
simulation_params.initialization = @(x) full_init_KdV(x);

[t_list,u_list] = PDE_solve(simulation_params);
clear('simulation_params')



%Markov solution

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 10;   %end of simulation
simulation_params.howoften = 100;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = N;          %number of positive modes to simulate

%full model with no approximations
simulation_params.initialization = @(x) full_init_KdV(x);

[t_markov,u_markov] = PDE_solve(simulation_params);
clear('simulation_params')




[x,exact] = make_real_space(u_list(1:N,1:length(t_list)),N);
[~,markov] = make_real_space(u_markov(1:N,1:length(t_list)),N);



%create animation comparing all three
figure
for i = 1:length(t_list)
    plot(x,exact(:,i),'linewidth',2)
    hold on
    plot(x,markov(:,i),'r','linewidth',2)
    axis([0,2*pi*(N-1)/N,-3,3])
    title(sprintf('t = %.1f', t_list(i)))
    legend('Exact','Markov model','location','southeast')
    drawnow
    hold off
    saveas(gcf,sprintf('markov_anim%i',i),'png')
end


[x,exact] = make_real_space(u_list,128);



%create animation comparing all three
figure
for i = 1:length(t_list)
    plot(x,exact(:,i),'linewidth',2)
    axis([0,2*pi*(128-1)/128,-3,3])
    title(sprintf('t = %.1f', t_list(i)))
    drawnow
    hold off
    saveas(gcf,sprintf('exact_anim%i',i),'png')
end