clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

N_full = 48;
alpha = 1;
N_list_coeff = 4:2:24;
tau_list = [0.0;0.2;0.4;0.6;0.8;1.0];
N_list = 6:2:14;

% Load calculated coefficients
% This data is generated using c_vs_N_tests_coeff_data.m
load(['c1_op_n_1_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c2_op_n_2_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c3_op_n_3_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])

% Choose tau value
tau = 1.0;

% Find coefficients corresponding to tau and N
ind = find(tau == tau_list);
c1_op_tau = c1_op{1,ind};
c2_op_tau = c2_op{1,ind};
c3_op_tau = c3_op{1,ind};
c4_op_tau = c4_op{1,ind};

ind_N_l = find(N_list_coeff == N_list(1));
ind_N_h = find(N_list_coeff == N_list(end));
c1_op_N = c1_op_tau(:,ind_N_l:ind_N_h);
c2_op_N = c2_op_tau(:,ind_N_l:ind_N_h);
c3_op_N = c3_op_tau(:,ind_N_l:ind_N_h);
c4_op_N = c4_op_tau(:,ind_N_l:ind_N_h);

[c1_laws,r1] = create_scaling_laws(N_list(1:end),c1_op_N(1,:));
[c2_laws,r2] = create_scaling_laws(N_list(1:end),c2_op_N(1:2,:));
[c3_laws,r3] = create_scaling_laws(N_list(1:end),c3_op_N(1:3,:));
[c4_laws,r4] = create_scaling_laws(N_list(1:end),c4_op_N(1:4,:));

% Save scaling laws
save(['scaling_laws_n_1_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_x_100' num2str(100*tau) '.mat'],'c1_laws')
save(['scaling_laws_n_2_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_x_100' num2str(100*tau) '.mat'],'c2_laws')
save(['scaling_laws_n_3_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_x_100' num2str(100*tau) '.mat'],'c3_laws')
save(['scaling_laws_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_x_100' num2str(100*tau) '.mat'],'c4_laws')


coeff_array = zeros(4,length(N_list),4);
coeff_array(1,:,1) = c1_op_N(1,:);
coeff_array(1:2,:,2) = c2_op_N(1:2,:);
coeff_array(1:3,:,3) = c3_op_N(1:3,:);
coeff_array(1:4,:,4) = c4_op_N(1:4,:);

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
saveas(gcf,sprintf('Euler_N_list_%i_to_%i_coeff_fits_c1234',N_list(1),N_list(end)),'eps')


figure(2)
legend([dots(2,2) dots(2,3) dots(2,4)],{'$n$ = 2','$n$ = 3','$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Euler_N_list_%i_to_%i_coeff_fits_c234',N_list(1),N_list(end)),'eps')


figure(3)
legend([dots(3,3) dots(3,4)],{'$n$ = 3','$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Euler_N_list_%i_to_%i_coeff_fits_c34',N_list(1),N_list(end)),'eps')

figure(4)
legend([dots(4,4)],{'$n$ = 4'},'location','northeast')
saveas(gcf,sprintf('Euler_N_list_%i_to_%i_coeff_fits_c4',N_list(1),N_list(end)),'eps')
