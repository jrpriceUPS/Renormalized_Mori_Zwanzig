clear all; close all

load slopes
load slopes2
load turn_times
load ens_max
load ens_max_time
load vort_max
load vort_max_time

N_list = 4:2:24;

figure(1)
hold off
plot(N_list,turn_times,'b.','markersize',20)
hold on
plot(N_list,ens_max_time,'r.','markersize',20)
plot(N_list,vort_max_time,'k.','markersize',20)
xlabel('N','fontsize',16)
ylabel('Time of events','fontsize',16)
legend('Energy decay begins','Maximum enstrophy attained','Maximum vorticity attained','location','southeast')
title('Event Times','fontsize',16)
saveas(gcf,'turn_times','png')

figure(2)
hold off
plot(N_list,slopes,'.','markersize',20)
xlabel('N','fontsize',16)
ylabel('Slope','fontsize',16)
title('Initial Decay Slope','fontsize',16)
saveas(gcf,'slopes','png')

figure(3)
hold off
plot(N_list,slopes2,'.','markersize',20)
xlabel('N','fontsize',16)
ylabel('Slope','fontsize',16)
title('Secondary Decay Slope','fontsize',16)
saveas(gcf,'slopes2','png')

figure(4)
hold off
plot(N_list,ens_max,'.','markersize',20)
xlabel('N','fontsize',16)
ylabel('Maximum enstrophy','fontsize',16)
title('Maximum Enstrophy','fontsize',16)
saveas(gcf,'enstrophy','png')

figure(5)
hold off
plot(N_list,vort_max,'.','markersize',20)
xlabel('N','fontsize',16)
ylabel('Maximum vorticity','fontsize',16)
title('Maximum Vorticity','fontsize',16)
saveas(gcf,'vorticity','png')


figure(6)
hold off
plot(N_list,ens_max,'.','markersize',20)
hold on
plot(N_list,vort_max,'r.','markersize',20)
xlabel('N','fontsize',16)
ylabel('Value','fontsize',16)
title('Maximum enstrophy and vorticity','fontsize',16)
legend('Maximum enstrophy','Maximum vorticity','location','northwest')
saveas(gcf,'maxima','png')