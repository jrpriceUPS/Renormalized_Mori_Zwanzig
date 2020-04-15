clear all;close all
N_list = 4:2:18;
colors = linspecer(length(N_list),'qualitative');
end_time = 1000;

for i = 1:length(N_list)
    
    N = N_list(i);
    full_legend{i} = sprintf('Fourth order N = %i ROM',N);
    
end

figure
hold on
for i = 1:length(N_list)
    
    for j = 1:i
        leg_sw{j} = full_legend{j};
    end
    
    leg_sw{i+1} = 'location';
    leg_sw{i+2} = 'southwest';
    
    N = N_list(i);
    load(sprintf('energy_%i_%i.mat',N,end_time))
    load(sprintf('t4_%i_%i',N,end_time))
    plot(log(t4),log(energy),'linewidth',1.5,'color',colors(i,:))
    
    
    legend(leg_sw{:})
    xlabel('log(time)','fontsize',20)
    ylabel('log(energy)','fontsize',20)
    axis([log(t4(2)),max(log(t4)),-10,0])
    ax = gca;
    ax.FontSize = 16;
    
end