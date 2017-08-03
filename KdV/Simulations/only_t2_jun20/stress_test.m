%compares t^2-model to Markov model

clear all;close all;

addpath ../../nonlinear
addpath ../../simulation_functions
addpath ../../analysis

epsilon = fliplr(0.065:0.005:0.1);
N_list = 8:2:40;

endtime = 100;
howoften = 100;

time_exact_array = zeros(length(N_list),length(epsilon));
time_markov_array = zeros(length(N_list),length(epsilon));
time_4_array = zeros(length(N_list),length(epsilon));
time_2_array = zeros(length(N_list),length(epsilon));

rel_err_markov_array = zeros(length(N_list),(endtime/1e-3)/howoften + 1,length(epsilon));
rel_err_4_array = zeros(length(N_list),(endtime/1e-3)/howoften + 1,length(epsilon));
rel_err_2_array = zeros(length(N_list),(endtime/1e-3)/howoften + 1,length(epsilon));

load coeff_array

for i = 1:length(epsilon)
    
    [timing_exact,timing_markov,timing_4,timing_2,rel_err_markov,rel_err_4,rel_err_2] = compare_sim_methods_just_t2(epsilon(i),coefficients,N_list,128,endtime,howoften);
    
    time_exact_array(:,i) = timing_exact;
    time_markov_array(:,i) = timing_markov;
    time_4_array(:,i) = timing_4;
    time_2_array(:,i) = timing_2;
    
    rel_err_markov_array(:,:,i) = rel_err_markov;
    rel_err_4_array(:,:,i) = rel_err_4;
    rel_err_2_array(:,:,i) = rel_err_2;
    
    epsilon(i)
    sum(time_exact_array(:,i)+time_markov_array(:,i)+time_4_array(:,i)+time_2_array(:,i))
    
end