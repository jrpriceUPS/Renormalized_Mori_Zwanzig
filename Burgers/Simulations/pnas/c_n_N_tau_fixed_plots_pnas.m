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

N_full = 16384;
alpha = 1;
epsilon = 1;
endtime = 1;
N = 14;
tau_list = 0.33;
t_slope_start = 15;
t_slope_end = 500;

% Create the "full" spectral solution data or load it
if ~(exist(sprintf(['u_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat']),'file') == 2)
    
    [t_list,u_list] = create_data_spec(alpha,N_full,endtime,epsilon);
    
    save(['u_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'u_list','-v7.3')
    save(['t_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'t_list','-v7.3')
    
else
    
    load(['u_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'u_list')
    load(['t_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'t_list')
    
end

% Find the resolved snapshots & trim the time and solution arrays or load
% the snapshot listing - this can be uncommented once we settle on final
% resolved criteria
if ~(exist(sprintf(['viable_snapshots_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat']),'file') == 2)
    
    [viable_snapshots,tmodel_size_list,tmodel_size_list_full] = resolve_array(u_list,t_list,alpha);
    
    save(['viable_snapshots_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'viable_snapshots');
    save(['tmodel_size_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'tmodel_size_list');
    save(['tmodel_size_list_full_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'tmodel_size_list_full');
    
else
    
    load(['viable_snapshots_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'viable_snapshots');
    load(['tmodel_size_list_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'tmodel_size_list');
    load(['tmodel_size_list_full_M_' num2str(N_full) '_endtime_' num2str(endtime) '_inveps_' num2str(1/epsilon) '.mat'],'tmodel_size_list_full');
    
end

% Trim the solution arrays
u_list = u_list(:,viable_snapshots);
t_list = t_list(viable_snapshots);

% Calculate the renromalization coefficients and corresponding error for
% tau_list
[c1_data,c2_data,c3_data,c4_data] = renormalize(alpha,N,u_list,t_list,tau_list);

% Extract variables from renormalize data structure
c1_op = c1_data.c1_op;
c2_op = c2_data.c2_op;
c3_op = c3_data.c3_op;
c4_op = c4_data.c4_op;

% Simulations and plots using tau corresponding to the minimum error from tau search
endtime = 1000;
howoften = 1000;
num_points = 10000;
dt = 1e-4;

% Create the exact upwind data or create it
if ~(exist(sprintf(['u_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat']),'file') == 2)
    [t_list,u_list,u_list_real] = create_data(alpha,num_points,endtime,dt,howoften,epsilon);
    
else
    load(['u_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'u_list');
    load(['t_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'t_list');
    load(['u_list_real_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'u_list_real');

end

% n = 1 ROM

simulation_params.N = N;
simulation_params.epsilon = epsilon;
simulation_params.alpha = alpha;
simulation_params.tau = c1_op(2,1);
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.time_range = t_list;
simulation_params.degree = 1;
simulation_params.coeffs = c1_op(1,1);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
[tc1B,uc1B] = ROM_PDE_solve(simulation_params);

energyc1B = get_energy(uc1B,N);

u_exact = u_list(1:N,:);

cfitROM1 = polyfit(log(tc1B(find(tc1B == t_slope_start):find(tc1B == t_slope_end))),log(energyc1B(find(tc1B == t_slope_start):find(tc1B == t_slope_end))), 1);

errc1B = ((get_energy(uc1B(:,1:length(tc1B)) - u_exact(:,1:length(tc1B)),N))./get_energy(u_exact(:,1:length(tc1B)),N)).';

even_log_space = exp(linspace(-2,log(tc1B(end)),50));
indexes = zeros(length(even_log_space),1);
for j = 1:length(even_log_space)
    [~,min_loc] = min(abs(tc1B - even_log_space(j)));
    indexes(j) = min_loc;
end

t1_e = tc1B(indexes);
energy1_e = energyc1B(indexes);
error1_e = errc1B(indexes);

% n = 2 ROM

clear u_exact even_log_space indexes

simulation_params.N = N;
simulation_params.epsilon = epsilon;
simulation_params.alpha = alpha;
simulation_params.tau = c2_op(3,1);
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.time_range = t_list;
simulation_params.degree = 2;
simulation_params.coeffs = c2_op(1:2,1);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
[tc2B,uc2B] = ROM_PDE_solve(simulation_params);

energyc2B = get_energy(uc2B,N);

u_exact = u_list(1:N,:);

cfitROM2 = polyfit(log(tc2B(find(tc2B == t_slope_start):find(tc2B == t_slope_end))),log(energyc2B(find(tc2B == t_slope_start):find(tc2B == t_slope_end))), 1);

errc2B = ((get_energy(uc2B(:,1:length(tc2B)) - u_exact(:,1:length(tc2B)),N))./get_energy(u_exact(:,1:length(tc2B)),N)).';

even_log_space = exp(linspace(-2,log(tc2B(end)),50));
indexes = zeros(length(even_log_space),1);
for j = 1:length(even_log_space)
    [~,min_loc] = min(abs(tc2B - even_log_space(j)));
    indexes(j) = min_loc;
end

t2_e = tc2B(indexes);
energy2_e = energyc2B(indexes);
error2_e = errc2B(indexes);

% n = 3 ROM

clear u_exact even_log_space indexes

simulation_params.N = N;
simulation_params.epsilon = epsilon;
simulation_params.alpha = alpha;
simulation_params.tau = c3_op(4,1);
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.time_range = t_list;
simulation_params.degree = 3;
simulation_params.coeffs = c3_op(1:3,1);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
[tc3B,uc3B] = ROM_PDE_solve(simulation_params);

energyc3B = get_energy(uc3B,N);

u_exact = u_list(1:N,:);

cfitROM3 = polyfit(log(tc3B(find(tc3B == t_slope_start):find(tc3B == t_slope_end))),log(energyc3B(find(tc3B == t_slope_start):find(tc3B == t_slope_end))), 1);

errc3B = ((get_energy(uc3B(:,1:length(tc3B)) - u_exact(:,1:length(tc3B)),N))./get_energy(u_exact(:,1:length(tc3B)),N)).';

even_log_space = exp(linspace(-2,log(tc3B(end)),50));
indexes = zeros(length(even_log_space),1);
for j = 1:length(even_log_space)
    [~,min_loc] = min(abs(tc3B - even_log_space(j)));
    indexes(j) = min_loc;
end

t3_e = tc3B(indexes);
energy3_e = energyc3B(indexes);
error3_e = errc3B(indexes);

% n = 4 ROM

clear u_exact even_log_space indexes

simulation_params.N = N;
simulation_params.epsilon = epsilon;
simulation_params.alpha = alpha;
simulation_params.tau = c4_op(5,1);
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.time_range = t_list;
simulation_params.degree = 4;
simulation_params.coeffs = c4_op(1:4,1);
simulation_params.initialization = @(x) ROM_init_Burgers(x);
[tc4B,uc4B] = ROM_PDE_solve(simulation_params);

energy_exact = get_energy(u_list,N);
energyc4B = get_energy(uc4B,N);

u_exact = u_list(1:N,:);

cfitexact = polyfit(log(t_list(find(t_list == t_slope_start):find(t_list == t_slope_end))).',log(energy_exact(find(t_list == t_slope_start):find(t_list == t_slope_end))), 1);
cfitROM4 = polyfit(log(tc4B(find(tc4B == t_slope_start):find(tc4B == t_slope_end))),log(energyc4B(find(tc4B == t_slope_start):find(tc4B == t_slope_end))), 1);

errc4B = ((get_energy(uc4B(:,1:length(tc4B)) - u_exact(:,1:length(tc4B)),N))./get_energy(u_exact(:,1:length(tc4B)),N)).';

even_log_space = exp(linspace(-2,log(tc4B(end)),50));
indexes = zeros(length(even_log_space),1);
for j = 1:length(even_log_space)
    [~,min_loc] = min(abs(tc4B - even_log_space(j)));
    indexes(j) = min_loc;
end

t4_e = tc4B(indexes);
energy4_e = energyc4B(indexes);
error4_e = errc4B(indexes);

leg{1} = ['First order $N$ = 14 ROM: slope = ' num2str(cfitROM1(1),'%3.2f')];
leg{2} = ['Second order $N$ = 14 ROM: slope = ' num2str(cfitROM2(1),'%3.2f')];
leg{3} = ['Third order $N$ = 14 ROM: slope = ' num2str(cfitROM3(1),'%3.2f')];
leg{4} = ['Fourth order $N$ = 14 ROM: slope = ' num2str(cfitROM4(1),'%3.2f')];

figure()
hold on
plot(log(t1_e),log(energy1_e),style{1},'DisplayName',leg{1})
plot(log(t2_e),log(energy2_e),style{2},'DisplayName',leg{2})
plot(log(t3_e),log(energy3_e),style{3},'DisplayName',leg{3})
plot(log(t4_e),log(energy4_e),style{4},'DisplayName',leg{4})
box on
xlim([-2,log(t4_e(end))])
ylim([-14,0])
xlabel('Log(Time)','fontsize',16)
ylabel('Log(Energy)','fontsize',16)

plot(log(t_list),log(energy_exact),'k--','linewidth',2,'DisplayName',['Second-order upwind solution: slope = ' num2str(cfitexact(1),'%3.2f')])

legend('location','southwest')
legend show
saveas(gcf,sprintf('Burgers_energy_multiple_orders_%i',N),'eps')

lege{1} = ['First order $N$ = 14 ROM'];
lege{2} = ['Second order $N$ = 14 ROM'];
lege{3} = ['Third order $N$ = 14 ROM'];
lege{4} = ['Fourth order $N$ = 14 ROM'];

figure()
hold on
plot(t1_e,error1_e,style{1},'DisplayName',lege{1})
plot(t2_e,error2_e,style{2},'DisplayName',lege{2})
plot(t3_e,error3_e,style{3},'DisplayName',lege{3})
plot(t4_e,error4_e,style{4},'DisplayName',lege{4})
box on
xlim([0,t1_e(end)])
xlabel('Time','fontsize',16)
ylabel('Error','fontsize',16)

legend('location','northwest')
legend show
saveas(gcf,sprintf('Burgers_error_multiple_orders_%i',N),'eps')