clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

[t_list,u_list] = upwind_burgers(1,10000,10,1e-4,100);

energy8 = get_energy(u_list,8);
figure
plot(t_list,energy8,'b','linewidth',2)
axis([0,10,0,0.55])
legend('Mass in first 8 Fourier modes','location','northeast')

saveas(gcf,'burger_energy','png')
close