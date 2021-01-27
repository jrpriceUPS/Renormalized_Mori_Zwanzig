clear all, close all, clc

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

style{1} = 'k-o';
style{2} = 'k-*';
style{3} = 'k-+';
style{4} = 'k-s';
style{5} = 'k-x';
style{6} = 'k-^';
style{7} = 'k->';
style{8} = 'k-<';
style{9} = 'k-p';
style{10} = 'k-h';

sim_endtime = 1000;
N_full = 48;
N = 14;
tau = 1;

leg{1} = ['First order $N$ = 14 ROM'];
leg{2} = ['Second order $N$ = 14 ROM'];
leg{3} = ['Third order $N$ = 14 ROM'];
leg{4} = ['Fourth order $N$ = 14 ROM'];

figure()
hold on
for i = 1:4
    
    clear u t t_e energy_e

    % This data is generated using c_vs_N_tests_sims.m
    if i == 4
    load(['t_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
    load(['u_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
    else
        load(['t_n_' num2str(i) '_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
        load(['u_n_' num2str(i) '_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
    end
    
    % Check that the sizes match
    if size(u,1) ~= N
        break
    end

    % Calculate the energy in the N modes of the ROM
    energy = get_3D_energy(u,N);
    
    even_log_space = exp(linspace(-2,log(t(end)),50));
    indexes = zeros(50,1);
    for j = 1:length(even_log_space)
        [~,min_loc] = min(abs(t - even_log_space(j)));
        indexes(j) = min_loc;
    end
    
    t_e = t(indexes);
    energy_e = energy(indexes);

    plot(log(t_e),log(energy_e),style{i},'DisplayName',leg{i})
    box on
    xlim([-2,log(t_e(end))])
    ylim([-11,0])
    xlabel('Log(Time)','fontsize',16)
    ylabel('Log(Energy)','fontsize',16)
    
end

hold off
legend('location','southwest')
legend show
saveas(gcf,sprintf('Euler_energy_multiple_orders_%i',N),'eps')
