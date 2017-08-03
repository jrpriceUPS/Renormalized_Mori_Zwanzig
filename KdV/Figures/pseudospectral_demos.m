%makes series of plots demonstrating pseudospectral methods
addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis


%exact solution
[t_exact,u_exact] = upwind_burgers(1,10000,10,0.0001,100);

%plot real space at different resolutions
[x,u] = make_real_space(u_exact,5000);
[x2,u2] = make_real_space(u_exact,8);
[x3,u3] = make_real_space(u_exact,16);

pseudospec = figure;
plot(x,u(:,175),'linewidth',2)
set(gca,'FontSize',16)
axis([0,2*pi,-1.25,1.25])
legend('Exact','location','northeast')
saveas(pseudospec,'exact_pseudospec','png')

hold on
plot(x2,u2(:,175),'r','linewidth',2)
axis([0,2*pi,-1.25,1.25])
legend('Exact','first 8 Fourier modes','location','northeast')
saveas(pseudospec,'pseudospec8','png')


plot(x3,u3(:,175),'k','linewidth',2)
axis([0,2*pi,-1.25,1.25])
legend('Exact','first 8 Fourier modes','first 16 Fourier modes','location','northeast')
saveas(pseudospec,'pseudospec16','png')
close

%show wigglyness when using many modes
[x4,u4] = make_real_space(u_exact,128);

wiggly_pseudo = figure;
plot(x4,u4(:,175),'linewidth',2)
set(gca,'FontSize',16)
axis([0,2*pi,-1.25,1.25])
saveas(wiggly_pseudo,'wiggly_pseudo','png')
close

%show energy flow for Burgers (monotonic)
burger_energy = figure;
plot(t_exact,get_energy(u_exact,8),'linewidth',2)
axis([0,10,0,0.55])
set(gca,'FontSize',16)
legend('Energy in first 8 Fourier modes')
saveas(burger_energy,'burger_energy','png')
close