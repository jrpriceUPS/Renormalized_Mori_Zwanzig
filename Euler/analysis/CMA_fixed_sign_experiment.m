clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

load u48
load u48_2

load t48
load t48_2

s = size(u);

u_both = zeros(s(1),s(2),s(3),s(4),s(5),length(t)+length(t2));
u_both(:,:,:,:,:,1:length(t)) = u;
u_both(:,:,:,:,:,length(t)+1:end) = u2;

t_both = [t;t2];

if ~(exist('tmodel_size_list48.mat','file') == 2)
    [tmodel_size_list,tmodel_size_list_full] = resolve_array(u_both,t_both);
    save('tmodel_size_list48','tmodel_size_list')
    save('tmodel_size_list_full48','tmodel_size_list_full')
end

load tmodel_size_list48
min_tol = 1e-16;
max_tol = 1e-10;

viable_snapshots = find(tmodel_size_list > min_tol & tmodel_size_list < max_tol);

% trim the arrays to those viable times
u_array = u_both(:,:,:,:,:,viable_snapshots);
t_array = t_both(viable_snapshots);

N_list = 4:2:24;

[dE,dE_CMA] = CMA_fixed_sign_check(u_array,N_list,t_array)