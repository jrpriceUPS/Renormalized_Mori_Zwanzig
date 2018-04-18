%clear all;close all

addpath ../../nonlinear
addpath ../../simulation_functions
addpath ../../analysis


N = 8;
alpha = 1;
dt = 1e-3;
endtime = 1.40492;
howoften = 10;

degree = 4;





simulation_params.N = N;
simulation_params.epsilon = 0;
simulation_params.alpha = alpha;
simulation_params.dt = dt;
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.howoften = 1;
simulation_params.tol = inf;
simulation_params.blowup = 1;
simulation_params.initial_condition = @(x) sin(x);


if degree == 0
    
    
    simulation_params_markov = simulation_params;
    simulation_params_markov.time_dependence = 1;
    simulation_params_markov.degree = 0;
    simulation_params_markov.initialization = @(x) ROM_init_Burgers(x);
    [t_Burgers,u_Burgers] = PDE_solve(simulation_params_markov);

    
    simulation_params_markov.name = 'full';
    simulation_params_markov.initialization = @(x) full_init_KdV(simulation_params_markov);
    [t_KdV,u_KdV] = PDE_solve2(simulation_params_markov);
    
end

if degree == 1
    
    simulation_params_tmodel = simulation_params;
    simulation_params_tmodel.time_dependence = 1;
    simulation_params_tmodel.degree = 1;
    simulation_params_tmodel.initialization = @(x) ROM_init_Burgers(x);
    [t_Burgers,u_Burgers] = PDE_solve(simulation_params_tmodel);
    
    simulation_params_tmodel.coeffs = 1;
    simulation_params_tmodel.initialization = @(x) complete_init_KdV(simulation_params_tmodel);
    [t_KdV,u_KdV] = PDE_solve2(simulation_params_tmodel);
    
end

if degree == 2
    simulation_params_t2model = simulation_params;
    simulation_params_t2model.time_dependence = 1;
    simulation_params_t2model.degree = 2;
    simulation_params_t2model.initialization = @(x) ROM_init_Burgers(x);
    [t_Burgers,u_Burgers] = PDE_solve(simulation_params_t2model);
    
    simulation_params_t2model.coeffs = [1;-1/2];
    simulation_params_t2model.initialization = @(x) complete_init_KdV(simulation_params_t2model);
    [t_KdV,u_KdV] = PDE_solve2(simulation_params_t2model);
    
    
end

if degree == 3
    
    simulation_params_t3model = simulation_params;
    simulation_params_t3model.time_dependence = 1;
    simulation_params_t3model.degree = 3;
    simulation_params_t3model.initialization = @(x) ROM_init_Burgers(x);
    [t_Burgers,u_Burgers] = PDE_solve(simulation_params_t3model);
    
    simulation_params_t3model.coeffs = [1;-1/2;1/6];
    simulation_params_t3model.initialization = @(x) complete_init_KdV(simulation_params_t3model);
    [t_KdV,u_KdV] = PDE_solve2(simulation_params_t3model);
    
end

if degree == 4
    
    simulation_params_t4model = simulation_params;
    simulation_params_t4model.time_dependence = 1;
    simulation_params_t4model.degree = 4;
    simulation_params_t4model.initialization = @(x) ROM_init_Burgers(x);
    [t_Burgers,u_Burgers] = PDE_solve(simulation_params_t4model);
    
    simulation_params_t4model.coeffs = [1;-1/2;1/6;-1/24];
    simulation_params_t4model.initialization = @(x) complete_init_KdV(simulation_params_t4model);
    [t_KdV,u_KdV] = PDE_solve2(simulation_params_t4model);
    
end





figure(2)
energy_Burgers = get_energy(u_Burgers,N);
energy_KdV = get_energy(u_KdV,N);
plot(log(t_Burgers),log(energy_Burgers),'b','linewidth',2)
hold on
plot(log(t_KdV),log(energy_KdV),'r','linewidth',2);
xlabel('log(time)','fontsize',16)
ylabel('log(energy)','fontsize',16)
legend('Burgers','KdV with epsilon = 0','location','southwest')