% produces plot to help visualize when the turn times converge to

clear all;close all;

load turn_times


plot(turn_times(1:end-1),diff(turn_times),'.','markersize',20)
ax = axis;
xlabel('T_N','fontsize',16)
ylabel('T_{N+2}-T_N','fontsize',16)
fit = polyfit(turn_times(1:end-1),diff(turn_times),1);
hold on
plot(ax(1:2),polyval(fit,ax(1:2)),'r');
plot([3.5,8],[0,0],'k')

intersect = -fit(2)/fit(1)

saveas(gcf,'turn_time_plot','png')