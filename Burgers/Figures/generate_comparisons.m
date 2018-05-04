function [times_array,energies_array,error_array] = generate_comparisons(N_list,endtime)
%
% [times_array,energies_array,error_array] = generate_comparisons(N_list,endtime)
%
% Compares the results of many Burgers' equation ROMs for many resolutions
% using calls to error_test
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%   N_list  =  set of resolutions
%
%  endtime  =  final time of simulations
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%     times_array  =  a cell array of times structures from error_test
%
%  energies_array  =  a cell array of energies structures from error_test
%
%     error_array  =  a cell array of errors structures from error_test

times_array = cell(length(N_list),1);
energies_array = cell(length(N_list),1);
error_array = cell(length(N_list),1);

for i = 1:length(N_list)
    [times,energies,errors] = error_test(N_list(i),1,endtime);
    times_array{i} = times;
    energies_array{i} = energies;
    error_array{i} = errors;
end