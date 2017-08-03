%runs simulation comparing Markov and order 4 ROM to exact solution

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

N_list = 20;
epsilon = 0.1;

%exact solution

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 100;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 128;          %number of positive modes to simulate

%full model with no approximations
simulation_params.name = 'full';

[t_list,u_list] = KdV_solve(simulation_params);
clear('simulation_params')



%Markov solution

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 100;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 20;          %number of positive modes to simulate

%full model with no approximations
simulation_params.name = 'full';

[t_markov,u_markov] = KdV_solve(simulation_params);
clear('simulation_params')


%ROM solution

simulation_params.epsilon = epsilon;  %coefficient on linear term in KdV
simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 100;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = inf;    %tolerance for identifying instabilities
simulation_params.N = 20;          %number of positive modes to simulate

%full model with no approximations
simulation_params.name = 'complete';

[t_ROM,u_ROM] = KdV_solve(simulation_params);


clear('simulation_params')

%plot energy evolution up to time 20

energy20 = figure;
set(gca,'FontSize',16)
hold on
plot(t_list(1:201),get_energy(u_list(:,1:201),20),'linewidth',2)
plot(t_markov(1:201),get_energy(u_markov(:,1:201),20),'r','linewidth',2)
plot(t_ROM(1:201),get_energy(u_ROM(:,1:201),20),'k','linewidth',2)
legend('Exact','Markov model','Complete ROM','location','northeast')
xlabel('time')
ylabel('Mass in resolved modes')
ax = axis;
axis([ax(1),ax(2),ax(3),0.5006])
saveas(energy20,'energy20','png')

[x,exact] = make_real_space(u_list(1:20,1:length(t_list)),20);
[~,markov] = make_real_space(u_markov(1:20,1:length(t_list)),20);
[~,ROM] = make_real_space(u_ROM(1:20,1:length(t_list)),20);

err_markov = real(sum((markov-exact).^2)./sum(exact.^2));
err_ROM = real(sum((ROM-exact).^2)./sum(exact.^2));


%plot real space error

error = figure;
hold on
set(gca,'FontSize',16)
plot(t_list,err_markov,'r','linewidth',2)
plot(t_list,err_ROM,'k','linewidth',2)
xlabel('time')
ylabel('global relative error')
legend('Markov model','Complete ROM','location','northwest')
saveas(error,'error','png')



%create animation comparing all three

for i = 1:5:length(t_list)
    plot(x,exact(:,i),'linewidth',2)
    hold on
    plot(x,markov(:,i),'r','linewidth',2)
    plot(x,ROM(:,i),'k','linewidth',2)
    axis([0,2*pi,-3,3])
    title(sprintf('t = %.1f', t_list(i)))
    legend('Exact','Markov model','Complete ROM','location','southwest')
    drawnow
    hold off
    saveas(gcf,sprintf('ROM_anim%i',(i-1)/5+1),'png')
end