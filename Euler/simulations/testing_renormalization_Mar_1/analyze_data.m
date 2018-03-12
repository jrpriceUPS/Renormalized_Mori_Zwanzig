function [coeff_array,scaling_laws] = analyze_data(N,min_tol,max_tol,time,print)
%
% A function to compute renormalization coefficients from a full simulation
% of size N using data where the flow out of the full model is between
% min_tol and max_tol using a least-squares fit
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%        N  =  resolution of full model  
%
%  min_tol  =  minimum energy flow out of t-model of full model (typically
%              0)
%
%  max_tol  =  maximum energy flow out of t-model of full model to use in
%              fit (typically 1e-16)
%
%     time  =  a logical variable (if 1, include time dependence, if 0 do
%              not include time dependence)
%
%    print  =  a logical variable (if 1, save plots, if 0 do not)
%
%
%%%%%%%%%
%OUTPUS:%
%%%%%%%%%
%
%   coeff_array  =  4 x (N - 2)/2 array of computed renormalization
%                   coefficients
%
%  scaling_laws  =  4 x 2 array of coefficients for fit scaling laws

% load relevant folders into the path
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

% load data
load(sprintf('t%i.mat',N))
load(sprintf('u%i.mat',N))
load(sprintf('tmodel_size_list%i.mat',N))

% isolate which timesteps fit the constraints
viable_snapshots = find(tmodel_size_list > min_tol & tmodel_size_list < max_tol);

% trim the arrays to those viable times
u_array = u(:,:,:,:,:,viable_snapshots);
t_array = t(viable_snapshots);

% compute the coefficients for all even modes smaller than the full
% simulation (except N = 2)
s = size(u);
N_list = 4:2:s(1)/2;

% compute the renormalization coefficients
[coeff_array,scaling_laws] = renormalize(u_array,N_list,t_array,time,print);