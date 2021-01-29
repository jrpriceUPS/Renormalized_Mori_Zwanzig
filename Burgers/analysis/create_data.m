function [t_list,u_list,u_list_real] = create_data(alpha,num_points,endtime,dt,howoften,epsilon)
%
%  create_data(alpha,num_points,endtime,dt,howoften)
%
%  Creates and saves exact solution data (u, t, u_real)
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       alpha  =  degree of nonlinearity
%
%  num_points  =  number of spatial points in upwind solution
%
%     endtime  =  final time of the simulation
%
%          dt  =  timestep
%
%    howoften  =  how often to save results (if 0, save all timesteps)
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%            t_list  =  times in the simulation
%
%            u_list  =  exact solution at each time
%
%       u_list_real  =  real space exact solution at each time

addpath ../../analysis
addpath ../../nonlinear
addpath ../../simulation_functions

% Compute and save exact solution
[t_list,u_list,u_list_real] = upwind_burgers(alpha,num_points,endtime,dt,howoften,epsilon);
save(['u_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'u_list','-v7.3');
save(['t_list_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'t_list','-v7.3');
save(['u_list_real_' num2str(endtime) '_num_points_' num2str(num_points) '_invdt_' num2str(1/dt) '_inveps_' num2str(1/epsilon) '.mat'],'u_list_real','-v7.3');