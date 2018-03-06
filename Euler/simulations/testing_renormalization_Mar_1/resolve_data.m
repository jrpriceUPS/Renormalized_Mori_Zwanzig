function resolve_data(N)

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

if exist(sprintf('tmodel_size_list%i.mat',N),'file') == 2
    
    load(sprintf('tmodel_size_list%i.mat',N))
    completed = length(tmodel_size_list);
    
else
    
    completed = 0;
    tmodel_size_list = [];
    
end

% load u_full_data
load(sprintf('u%i.mat',N))
load(sprintf('t%i.mat',N))

if length(t) > completed
    
    u_new = u(:,:,:,:,:,completed+1:end);
    t_new = t(completed+1:end);
    
    tmodel_size_list_new = resolve_array(u_new,t_new);
    
    tmodel_size_list = [tmodel_size_list;tmodel_size_list_new];
    save(sprintf('tmodel_size_list%i.mat',N),'tmodel_size_list')
    
end