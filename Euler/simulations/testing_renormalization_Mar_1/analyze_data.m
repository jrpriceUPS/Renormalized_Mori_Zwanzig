function [coeff_array,scaling_laws] = analyze_data(N,min_tol,max_tol,time,print)

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

load(sprintf('t%i.mat',N))
load(sprintf('u%i.mat',N))
load(sprintf('tmodel_size_list%i.mat',N))

viable_snapshots = find(tmodel_size_list > min_tol & tmodel_size_list < max_tol);

u_array = u(:,:,:,:,:,viable_snapshots);
t_array = t(viable_snapshots);

s = size(u);
N_list = 4:2:s(1)-2;

[coeff_array,scaling_laws] = renormalize(u_array,N_list,t_array,time,print);