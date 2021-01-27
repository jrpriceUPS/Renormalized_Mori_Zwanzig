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
style{10} = 'k-d';
style{11} = 'k-h';

N_full = 16384;
alpha = 1;
epsilon = 1;
endtime = 1;
N_list = 6:2:14;
tau_list = 0.40;
tau_fit = tau_list;
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
[c1_data,c2_data,c3_data,c4_data] = renormalize(alpha,N_list,u_list,t_list,tau_list);

% Extract variables from renormalize data structure
c1_op = c1_data.c1_op;
c2_op = c2_data.c2_op;
c3_op = c3_data.c3_op;
c4_op = c4_data.c4_op;

% Create the scaling laws for renormalization coefficients only
[c1_laws,r1] = create_scaling_laws(N_list,c1_op(1,:));
[c2_laws,r2] = create_scaling_laws(N_list,c2_op(1:2,:));
[c3_laws,r3] = create_scaling_laws(N_list,c3_op(1:3,:));
[c4_laws,r4] = create_scaling_laws(N_list,c4_op(1:4,:));

% % Save scaling laws
% save(['scaling_laws_n_1_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_inveps_' num2str(1/epsilon) '.mat'],'c1_laws')
% save(['scaling_laws_n_2_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_inveps_' num2str(1/epsilon) '.mat'],'c2_laws')
% save(['scaling_laws_n_3_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_inveps_' num2str(1/epsilon) '.mat'],'c3_laws')
% save(['scaling_laws_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_inveps_' num2str(1/epsilon) '.mat'],'c4_laws')
 
% load data into arrays for looping
coeff_array = zeros(4,length(N_list),4);
coeff_array(1,:,1) = c1_op(1,:);
coeff_array(1:2,:,2) = c2_op(1:2,:);
coeff_array(1:3,:,3) = c3_op(1:3,:);
coeff_array(1:4,:,4) = c4_op(1:4,:);

tau_array = zeros(1,length(N_list),4);
tau_array(1,:,1) = c1_op(2,:);
tau_array(1,:,2) = c2_op(3,:);
tau_array(1,:,3) = c3_op(4,:);
tau_array(1,:,4) = c4_op(5,:);

scaling_laws = zeros(4,2,4);
scaling_laws(1,:,1) = c1_laws;
scaling_laws(1:2,:,2) = c2_laws;
scaling_laws(1:3,:,3) = c3_laws;
scaling_laws(1:4,:,4) = c4_laws;

% Plot the scaling laws
style{1} = 'k*';
size(1) = 5;
style{2} = 'ko';
size(2) = 5;
style{3} = 'ks';
size(3) = 5;
style{4} = 'k.';
size(4) = 20;
dots = zeros(4,4);
colors = hsv(4);

% plot the data and fits
for j = 1:4
    figure(1)
    a = plot(log(N_list),log(squeeze(coeff_array(1,:,j))),style{j},'markersize',size(j));
    dots(1,j) = a;
    hold on
    plot([log(N_list(1))*0.95,log(N_list(end))*1.05],polyval(scaling_laws(1,:,j),[log(N_list(1))*0.95,log(N_list(end))*1.05]),'color','k')
    xlim([log(N_list(1))*0.95,log(N_list(end))*1.05])
    xlabel('Log($N$)','fontsize',16)
    ylabel('Log($a_1$)','fontsize',16)
    
    if j > 1
        
        figure(2)
        a = plot(log(N_list),log(-squeeze(coeff_array(2,:,j))),style{j},'markersize',size(j));
        dots(2,j) = a;
        hold on
        plot([log(N_list(1))*0.95,log(N_list(end))*1.05],polyval(scaling_laws(2,:,j),[log(N_list(1))*0.95,log(N_list(end))*1.05]),'color','k')
        xlim([log(N_list(1))*0.95,log(N_list(end))*1.05])
        xlabel('Log($N$)','fontsize',16)
        ylabel('Log($-a_2$)','fontsize',16)
        
        if j > 2
            
            figure(3)
            a = plot(log(N_list),log(squeeze(coeff_array(3,:,j))),style{j},'markersize',size(j));
            dots(3,j) = a;
            hold on
            plot([log(N_list(1))*0.95,log(N_list(end))*1.05],polyval(scaling_laws(3,:,j),[log(N_list(1))*0.95,log(N_list(end))*1.05]),'color','k')
            xlim([log(N_list(1))*0.95,log(N_list(end))*1.05])
            xlabel('Log($N$)','fontsize',16)
            ylabel('Log($a_3$)','fontsize',16)
            
            if j > 3
                
                figure(4)
                a = plot(log(N_list),log(-squeeze(coeff_array(4,:,j))),style{j},'markersize',size(j));
                dots(4,j) = a;
                hold on
                plot([log(N_list(1))*0.95,log(N_list(end))*1.05],polyval(scaling_laws(4,:,j),[log(N_list(1))*0.95,log(N_list(end))*1.05]),'color','k')
                xlim([log(N_list(1))*0.95,log(N_list(end))*1.05])
                xlabel('Log($N$)','fontsize',16)
                ylabel('Log($-a_4$)','fontsize',16)
            end
        end
    end
end

figure(1)
legend([dots(1,1) dots(1,2) dots(1,3) dots(1,4)],{'$n$ = 1','$n$ = 2','$n$ = 3','$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Burgers_N_list_%i_to_%i_coeff_fits_c1234',N_list(1),N_list(end)),'eps')


figure(2)
legend([dots(2,2) dots(2,3) dots(2,4)],{'$n$ = 2','$n$ = 3','$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Burgers_N_list_%i_to_%i_coeff_fits_c234',N_list(1),N_list(end)),'eps')


figure(3)
legend([dots(3,3) dots(3,4)],{'$n$ = 3','$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Burgers_N_list_%i_to_%i_coeff_fits_c34',N_list(1),N_list(end)),'eps')

figure(4)
legend([dots(4,4)],{'$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Burgers_N_list_%i_to_%i_coeff_fits_c4',N_list(1),N_list(end)),'eps')

% Simulations and plots using tau corresponding to the minimum error from tau search
% endtime = 1000;
endtime = 1000;
howoften = 1000;
num_points = 10000;
dt = 1e-4;
N_list = 18:4:34;

% Create the exact upwind data or create it
if ~(exist(sprintf(['u_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat']),'file') == 2)
    [t_list,u_list,u_list_real] = create_data(alpha,num_points,endtime,dt,howoften,epsilon);
    
else
    load(['u_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'u_list');
    load(['t_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'t_list');
    load(['u_list_real_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'u_list_real');

end

% Create empty dictionaries
tc4B{1,length(N_list)} = [];
uc4B{1,length(N_list)} = [];
energyc4B{1,length(N_list)} = [];
errc4B{1,length(N_list)} = [];
R0{1,length(N_list)} = [];
R1{1,length(N_list)} = [];
R2{1,length(N_list)} = [];
R3{1,length(N_list)} = [];
R4{1,length(N_list)} = [];
RTS{1,length(N_list)} = [];

% Preallocate arrays
errc4B_val = zeros(length(N_list),1);

% n = 4 ROMs
for i = 1:length(N_list)
    
    N = N_list(i);
    coeffs = zeros(5,1);
    
    coeffs(1) = exp(c4_laws(1,2))*N^c4_laws(1,1);
    coeffs(2) = -exp(c4_laws(2,2))*N^c4_laws(2,1);
    coeffs(3) = exp(c4_laws(3,2))*N^c4_laws(3,1);
    coeffs(4) = -exp(c4_laws(4,2))*N^c4_laws(4,1);
    coeffs(5) = tau_fit;
    
    simulation_params.N = N;
    simulation_params.epsilon = epsilon;
    simulation_params.alpha = alpha;
    simulation_params.tau = coeffs(5);
    simulation_params.endtime = endtime;
    simulation_params.print_time = 1;
    simulation_params.time_range = t_list;
    simulation_params.degree = 4;
    simulation_params.coeffs = coeffs(1:4);
    simulation_params.initialization = @(x) ROM_init_Burgers(x);
    [tc4B{i},uc4B{i}] = ROM_PDE_solve(simulation_params);
    
    energy_exact = get_energy(u_list,N);
    energyc4B{i} = get_energy(uc4B{i},N);
    
    u_exact = u_list(1:N,:);
    errc4B{i} = ((get_energy(uc4B{i}(:,1:length(tc4B{i})) - u_exact(:,1:length(tc4B{i})),N))./get_energy(u_exact(:,1:length(tc4B{i})),N)).';
    
    
    if length(tc4B{i}) > 0
        errc4B_val(i) = (1/tc4B{i}(end))*trapz(tc4B{i},errc4B{i});
    else
        errc4B_val(i) = 0;
    end
    
end

cfitexact = polyfit(log(t_list(find(t_list == t_slope_start):find(t_list == t_slope_end))).',log(energy_exact(find(t_list == t_slope_start):find(t_list == t_slope_end))), 1);

figure()
hold on
for i = 1:length(N_list)
    
    cfitROM = polyfit(log(tc4B{i}(find(tc4B{i} == t_slope_start):find(tc4B{i} == t_slope_end))),log(energyc4B{i}(find(tc4B{i} == t_slope_start):find(tc4B{i} == t_slope_end))), 1);
    
    even_log_space = exp(linspace(-2,log(tc4B{i}(end)),50));
    indexes = zeros(length(even_log_space),1);
    for j = 1:length(even_log_space)
        [~,min_loc] = min(abs(tc4B{i} - even_log_space(j)));
        indexes(j) = min_loc;
    end
    
    t_e = tc4B{i}(indexes);
    energy_e = energyc4B{i}(indexes);
    
    txt = ['Fourth order N = ' num2str(N_list(i)) ' ROM: slope = ' num2str(cfitROM(1),'%3.2f')];
    plot(log(t_e),log(energy_e),style{5+i},'DisplayName',txt)
    box on
    xlim([-2,log(t_e(end))])
    ylim([-14,0])
    xlabel('Log(Time)','fontsize',16)
    ylabel('Log(Energy)','fontsize',16)
    
    clear t_e energy_e
    
end

plot(log(t_list),log(energy_exact),'k--','linewidth',2,'DisplayName',['Second-order upwind solution: slope = ' num2str(cfitexact(1),'%3.2f')])

hold off
legend('location','southwest')
saveas(gcf,sprintf('Burgers_energy_scaling_laws_multiple_N_%i_to_%i',N_list(1),N_list(end)),'eps')

figure()
hold on
for i = 1:length(N_list)
    
    even_space = linspace(0,tc4B{i}(end),50);
    indexes = zeros(length(even_space),1);
    for j = 1:length(even_space)
        [~,min_loc] = min(abs(tc4B{i} - even_space(j)));
        indexes(j) = min_loc;
    end
    
    t_e = tc4B{i}(indexes);
    error_e = errc4B{i}(indexes);

    txt = ['Fourth order N = ' num2str(N_list(i)) ' ROM'];
    plot(t_e,error_e,style{5+i},'DisplayName',txt)
    box on
    xlim([0,t_e(end)])
    ylim([0,0.2])
    xlabel('Time','fontsize',16)
    ylabel('Error','fontsize',16)
    
    clear t_e error_e
        
end

hold off
legend('location','northwest')
saveas(gcf,sprintf('Burgers_error_scaling_laws_multiple_N_%i_to_%i',N_list(1),N_list(end)),'eps')