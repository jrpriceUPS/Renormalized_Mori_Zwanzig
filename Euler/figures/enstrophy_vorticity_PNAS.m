clear all;close all;

load ens_max
load vort_max

N_list = 4:2:24;

plot(N_list,ens_max,'k.','markersize',20)
hold on
plot(N_list,vort_max,'k*','markersize',10)
xlabel('N','fontsize',16)
ylabel('Maximum value','fontsize',16)
title('Maximum enstrophy and vorticity','fontsize',16)
legend('Enstrophy','Vorticity','location','northwest')

saveas(gcf,'maxima','png')
close all