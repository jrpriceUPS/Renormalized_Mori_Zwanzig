function compute_vorticities(N_list)

format long
close all

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

colors = linspecer(length(N_list));

figure(1)
hold on
figure(2)
hold on


for i = 1:length(N_list)
    
    N = N_list(i);
    full_legend{i} = sprintf('Fourth order N = %i ROM',N);
    
    load(sprintf('u_array4_%i.mat',N))
    load(sprintf('t4_%i',N))
    
    [v1,v2] = vorticity(u_array4);
    
    figure(1)
    plot(t4,v1,'linewidth',2,'color',colors(i,:))
    
    figure(2)
    plot(t4,v2,'linewidth',2,'color',colors(i,:))
    
end

figure(1)
title('Infinity norm of vorticity','fontsize',16)
legend(full_legend)
xlabel('time','fontsize',16)
ylabel('Maximum component of vorticity','fontsize',16)
saveas(gcf,'vorticity1','png')

figure(2)
title('Maximal 2-norm of vorticity','fontsize',16)
legend(full_legend)
xlabel('time','fontsize',16)
ylabel('Maximum 2-norm of vorticity','fontsize',16)
saveas(gcf,'vorticity2','png')
