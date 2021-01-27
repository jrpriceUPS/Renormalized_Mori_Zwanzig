clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

N_full = 48;
end_time = 1.5;
alpha = 1;
N_list = 4:2:24;
tau_list = [0.0;0.2;0.4;0.6;0.8;1.0];

c1_op{1,length(tau_list)} = [];
E1_op{1,length(tau_list)} = [];
Axcond1_op{1,length(tau_list)} = [];

c2_op{1,length(tau_list)} = [];
E2_op{1,length(tau_list)} = [];
Axcond2_op{1,length(tau_list)} = [];

c3_op{1,length(tau_list)} = [];
E3_op{1,length(tau_list)} = [];
Axcond3_op{1,length(tau_list)} = [];

c4_op{1,length(tau_list)} = [];
E4_op{1,length(tau_list)} = [];
Axcond4_op{1,length(tau_list)} = [];

for i = 1:length(tau_list)

    % Create the "full" spectral solution data or load it
    if ~(exist(sprintf(['u_' num2str(N_full) '_endtime_x10_' num2str(end_time*10) '.mat']),'file') == 2)
        
        [u,t] = create_data(N_full,end_time,alpha);
            
    else
        
        load(['u_' num2str(N_full) '_endtime_x10_' num2str(end_time*10) '.mat'],'u')
        load(['t_' num2str(N_full) '_endtime_x10_' num2str(end_time*10) '.mat'],'t')
        
    end
    
    % Calculate the energy flowing out the t-model for N_full/2 and N_full
    if ~(exist(sprintf(['tmodel_size_list_' num2str(N_full) '.mat']),'file') == 2)
        
        [tmodel_size_list,tmodel_size_list_full] = resolve_array(u,t);
        
        save(['tmodel_size_list_' num2str(N_full)],'tmodel_size_list')
        save(['tmodel_size_list_full_' num2str(N_full)],'tmodel_size_list_full')
        
    else
        
        load(['tmodel_size_list_' num2str(N_full)],'tmodel_size_list')
        load(['tmodel_size_list_full_' num2str(N_full)],'tmodel_size_list_full')
        
    end
    
    % Find the viable snapshots
    max_tol = 1e-10;
    viable_snapshots = find(tmodel_size_list < max_tol);

    % Trim the arrays to those viable times
    u_array = u(:,:,:,:,:,viable_snapshots);
    t_array = t(viable_snapshots);
    
    [c1_data,c2_data,c3_data,c4_data] = renormalize(u_array,N_list,t_array,tau_list(i));
    
    c1_op{i} = c1_data.c1_op;
    E1_op{i} = c1_data.E1_op;
    Axcond1_op{i} = c1_data.Axcond1_op;
    
    c2_op{i} = c2_data.c2_op;
    E2_op{i} = c2_data.E2_op;
    Axcond2_op{i} = c2_data.Axcond2_op;
    
    c3_op{i} = c3_data.c3_op;
    E3_op{i} = c3_data.E3_op;
    Axcond3_op{i} = c3_data.Axcond3_op;
    
    c4_op{i} = c4_data.c4_op;
    E4_op{i} = c4_data.E4_op;
    Axcond4_op{i} = c4_data.Axcond4_op;

end

save(['c1_op_n_1_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'],'c1_op')
save(['E1_op_n_1_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'E1_op')
save(['Axcond1_op_n_1_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'Axcond1_op')

save(['c2_op_n_2_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'],'c2_op')
save(['E2_op_n_2_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'E2_op')
save(['Axcond2_op_n_2_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'Axcond2_op')

save(['c3_op_n_3_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'],'c3_op')
save(['E3_op_n_3_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'E3_op')
save(['Axcond3_op_n_3_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'Axcond3_op')

save(['c4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'],'c4_op')
save(['E4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'E4_op')
save(['Axcond4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end))  '.mat'],'Axcond4_op')