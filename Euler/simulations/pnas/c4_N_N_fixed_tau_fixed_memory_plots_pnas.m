clear all; close all; clc;

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
tau = 1.0;

% Load t and R data for specified N and tau
% This data is generated using c_vs_N_test_memory_plots.m
load(['t_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
load(['R0_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
load(['R1_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
load(['R2_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
load(['R3_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
load(['R4_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])
load(['RT_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'])

RT = (R1_test+R2_test+R3_test+R4_test).';
RTS = trapz(t,RT,1);

figure()

even_log_space = exp(linspace(-2,log(t(end)),150));
indexes = zeros(length(even_log_space),1);
for j = 1:length(even_log_space)
    [~,min_loc] = min(abs(t - even_log_space(j)));
    indexes(j) = min_loc;
end

t_e = t(indexes);
R1_e = R1_test(indexes);
R2_e = R2_test(indexes);
R3_e = R3_test(indexes);
R4_e = R4_test(indexes);

hold on
plot(log(t_e),R1_e,style{2},'DisplayName','$\frac{1}{2} \sum_{k \in F} \alpha_1(t) t^1 \Delta E^1_k(t)$')
plot(log(t_e),R2_e,style{3},'DisplayName','$\frac{1}{2} \sum_{k \in F} \alpha_2(t) t^2 \Delta E^2_k(t)$')
plot(log(t_e),R3_e,style{4},'DisplayName','$\frac{1}{2} \sum_{k \in F} \alpha_3(t) t^3 \Delta E^3_k(t)$')
plot(log(t_e),R4_e,style{5},'DisplayName','$\frac{1}{2} \sum_{k \in F} \alpha_4(t) t^4 \Delta E^4_k(t)$')
plot(log(t),RT,'k-.','linewidth',2,'DisplayName','$\frac{1}{2}  \sum_{i=1}^4 \sum_{k \in F} \alpha_i(t) t^i \Delta E^i_k(t)$')
box on
xlim([-2,log(t_e(end))])
ylim([-0.025,0.010])
xlabel('Log(Time)','fontsize',16)
ylabel('Contribution to energy rate of change','fontsize',16)
hold off
legend('location','southeast')


saveas(gcf,sprintf('Euler_memory_energy_N_%i_tau_%i',N,tau),'eps')

    

