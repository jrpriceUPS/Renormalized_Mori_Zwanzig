function [t_list,u_list] = create_data_spec(alpha,N,endtime,epsilon)
%
%  create_data_spec(alpha,N,endtime,epsilon)
%
%  Creates "full" solution data (u, t)
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%             alpha  =  degree of nonlinearity
%
%                 N  =  number of positive spectral modes
%
%           endtime  =  final time of the simulation
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%            t_list  =  times in the simulation
%
%            u_list  =  exact solution at each time

addpath ../../analysis
addpath ../../nonlinear
addpath ../../simulation_functions

simulation_params.epsilon = epsilon;
simulation_params.N = N;
M = 3*N;
simulation_params.M = M;
simulation_params.alpha = alpha;
simulation_params.endtime = endtime;
simulation_params.print_time = 1;
simulation_params.initial_condition = @(x) epsilon*sin(x);
simulation_params.initialization = @(x) full_init_Burgers(x);
[t_list,u_list] = full_PDE_solve(simulation_params);

t_list = t_list.';