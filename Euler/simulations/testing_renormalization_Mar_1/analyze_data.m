%clear all; close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

tol = inf;
N_list = [2,4,6,8];

% load u_full_data
u_full_data = u_array_full;

%tol = inf
%u = resolve_array(u,tol);


coeff_array = renormalize(u_full_data,N_list);