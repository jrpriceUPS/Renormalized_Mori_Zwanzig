close all

addpath ../../nonlinear
addpath ../../simulation_functions
addpath ../../analysis

N = 8;
alpha = 1;
dt = 1e-4;
num_points = 10000;
endtime = 10;
howoften = 10;

degree = 0;





simulation_params.N = N;
simulation_params.alpha = alpha;
simulation_params.dt = dt;
simulation_params.endtime = endtime;
simulation_params.print_time = 1;

%[t_full,u_full] = upwind_burgers(alpha,num_points,endtime,dt,howoften);

simulation_params_markov = simulation_params;
simulation_params_markov.time_dependence = 1;
simulation_params_markov.degree = 0;
simulation_params_markov.initialization = @(x) ROM_init_Burgers(x);
[t0,u0] = PDE_solve(simulation_params_markov);

simulation_params_tmodel = simulation_params;
simulation_params_tmodel.time_dependence = 1;
simulation_params_tmodel.degree = 1;
simulation_params_tmodel.initialization = @(x) ROM_init_Burgers(x);
[t1,u1] = PDE_solve(simulation_params_tmodel);

simulation_params_t2model = simulation_params;
simulation_params_t2model.time_dependence = 1;
simulation_params_t2model.degree = 2;
simulation_params_t2model.initialization = @(x) ROM_init_Burgers(x);
[t2,u2] = PDE_solve(simulation_params_t2model);

simulation_params_t3model = simulation_params;
simulation_params_t3model.time_dependence = 1;
simulation_params_t3model.degree = 3;
simulation_params_t3model.initialization = @(x) ROM_init_Burgers(x);
[t3,u3] = PDE_solve(simulation_params_t3model);

simulation_params_t4model = simulation_params;
simulation_params_t4model.time_dependence = 1;
simulation_params_t4model.degree = 4;
simulation_params_t4model.initialization = @(x) ROM_init_Burgers(x);
[t4,u4] = PDE_solve(simulation_params_t4model);






energy_full = get_energy(u_full,N);
energy0 = get_energy(u0,N);
energy1 = get_energy(u1,N);
energy2 = get_energy(u2,N);
energy3 = get_energy(u3,N);
energy4 = get_energy(u4,N);
plot(log(t_full),log(energy_full),'b','linewidth',2)
hold on
plot(log(t0),log(energy0),'r','linewidth',2);
plot(log(t1),log(energy1),'k','linewidth',2);
plot(log(t2),log(energy2),'g','linewidth',2);
plot(log(t3),log(energy3),'m','linewidth',2);
plot(log(t4),log(energy4),'c','linewidth',2);
axis([log(t_full(2)),log(t_full(end)),-3,0])
title(sprintf('Energy in resolved modes N = %i, Burgers',N),'fontsize',16)
xlabel('log(time)','fontsize',16)
ylabel('log(energy)','fontsize',16)
legend('Exact (upwind method)','Markov model','t-model','t^2-model','t^3-model','t^4-model','location','southwest')