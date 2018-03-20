function resolve_data(N)
%
% A function to compute the hypothetical flow of energy out of the full
% model of size N
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%         N  =  resolution of full model
%
%
%%%%%%%%%
%OUTPUS:%
%%%%%%%%%
%
%  tmodel_size_listN.dat  =  saved vector of the time derivative of the
%                            t-model of the full model at associated times

% load relevant folders into the path
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

% load data
load(sprintf('u%i.mat',N))
load(sprintf('t%i.mat',N))

% use the resolve_array function to resolve each term
[tmodel_size_list,tmodel_size_list_full] = resolve_array(u,t);

% append the results and save it to the directory
save(sprintf('tmodel_size_list%i.mat',N),'tmodel_size_list')
save(sprintf('tmodel_size_list_full%i.mat',N),'tmodel_size_list_full')