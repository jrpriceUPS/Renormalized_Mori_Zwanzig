clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

n_list = [12:1:17];
N_full_list = 2.^n_list;
alpha = 1;
epsilon = 1;
endtime = 1;
N_list = 6:2:14;
tau_list = [-1:0.01:1];

for i = 1:length(N_full_list)
    
    N_full = N_full_list(i);
    
    [t_list,u_list] = create_data_spec(alpha,N_full,endtime,epsilon);
    
    [viable_snapshots,tmodel_size_list,tmodel_size_list_full] = resolve_array(u_list,t_list,alpha);
    
    % Trim the solution arrays
    u_list = u_list(:,viable_snapshots);
    t_list = t_list(viable_snapshots);
    
    % Calculate the renromalization coefficients and corresponding error for
    % tau_list
    [c1_data,c2_data,c3_data,c4_data] = renormalize(alpha,N_list,u_list,t_list,tau_list);
    
    % Extract variables from renormalize data structure
    c4_loc = c4_data.c4_op;
    
    % Save coefficient array
    save(['coeff_array_n_4_M_' num2str(N_full_list(i)) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_inveps_' num2str(1/epsilon) '_dt_5e5_.mat'],'c4_loc')

end


