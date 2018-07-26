addpath ../analysis

increment = 0.001;
max_percent = 0.995;

percentages = increment:increment:max_percent;
chart = zeros(length(N_list),length(percentages));

for i = 1:length(percentages)
    chart(:,i) = decay_begins(N_list,percentages(i),1000);
end

colors = linspecer(length(N_list),'qualitative');

leg = cell(length(N_list),1);

figure(1)
hold off
for i = 1:length(N_list)
    plot(chart(i,:),percentages,'linewidth',1.5,'color',colors(i,:))
    hold on;
    leg{i} = sprintf('N = %i fourth order ROM',N_list(i));
end

xlabel('Time','fontsize',16)
ylabel('Percentage drained','fontsize',16)
legend(leg{:},'location','southeast');
axis([0,50,0,max_percent])
saveas(gcf,'percentages','png')
