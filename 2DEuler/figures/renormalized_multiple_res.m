function [slopes,slopes2,turn_times,ens_max,ens_max_time,vort_max,vort_max_time] = renormalized_multiple_res(N_list,end_time)
%
%  [slopes,slopes2,turn_times,ens_max,ens_max_time,vort_max,vort_max_time] = renormalized_multiple_res(N_list,end_time,filetype)
%
%  A function to solve Euler's equations with 4th order complete memory
%  approximation ROMs for a variety of resolutions and plot the results
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%         N  =  resolution
%
%  end_time  =  time to run simulation to
%
%  filetype  =  file type in which to save figures
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%         slopes  =  the initial log-log energy decay slopes
% 
%        slopes2  =  the secondary log-log energy decay slopes
%
%     turn_times  =  the times at which energy begins draining from the system
%
%        ens_max  =  the maximum enstrophy
%
%   ens_max_time  =  the time at which the maximum enstrophy is achieved
%
%       vort_max  =  the maximum of the vorticity
%
%  vort_max_time  =  the time at which the maximum vorticity is achieved

%initialize everything needed for the simulation
format long
close all

% addpath ../simulation_functions
% addpath ../nonlinear
% addpath ../analysis

%colors = linspecer(length(N_list),'qualitative');

slopes = zeros(length(N_list),1);
slopes2 = zeros(length(N_list),1);
turn_times = zeros(length(N_list),1);
ens_max = zeros(length(N_list),1);
ens_max_time = zeros(length(N_list),1);
vort_max = zeros(length(N_list),1);
vort_max_time = zeros(length(N_list),1);

load('u64.mat')
load('t64.mat')

u_full = zeros(128,128,2,length(t));
for i = 1:length(t)
    u_full(:,:,:,i) = u_fullify(u(:,:,:,:,i),64);
end

for N = N_list
    
    M = 3*N;
    
    % uniform grid
    x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
    y = x;
    
    % get initial condition
    u = zeros(N,N,2,2,length(t));
    for i = 1:length(t)
        u(:,:,:,:,i) = u_squishify(u_full(:,:,:,i),N);
    end
    u_init = u(:,:,:,:,1);
    init_energy = get_2D_energy(u_init,N);
    
    % make k array
    k_vec = [0:M-1,-M:1:-1];
    [kx,ky] = ndgrid(k_vec,k_vec);
    k = zeros(2*M,2*M,2);
    k(:,:,1) = kx;
    k(:,:,2) = ky;
    
    params.k = k;
    params.N = N;
    params.M = M;
    params.a = 2:M;
    params.b = 2*M:-1:M+2;
    params.a_tilde = N+1:M;
    params.a_tilde2 = 2*N+1:M;
    params.print_time = 1;
    params.no_time = 0;
    
    detect_blowup = @(T,Y) energy_grow(T,Y,init_energy,1,params);
    
    params4 = params;
    params4.func = @(x) t4model_RHS(x);
    params4.coeff = [1 -1/2 1/6 -1/24];
        
        % run the simulation
        options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3,'Events',detect_blowup);
        [t4,u_raw4] = ode45(@(t,u) RHS(u,t,params4),[0,end_time],u_init(:),options);
        
        
        % reshape the output array into an intelligible shape (should make this a
        % separate function later)
        u_array4 = zeros([size(u_init) length(t4)]);
        for i = 1:length(t4)
            u_array4(:,:,:,:,i) = reshape(u_raw4(i,:),[N,N,2,2]);
        end
        
        save(sprintf('t4_%i_%i',N,end_time),'t4');
        save(sprintf('u_array4_%i_%i.mat',N,end_time),'u_array4');
    
    params3 = params;
    params3.func = @(x) t3model_RHS(x);
    params3.coeff = [1 -1/2 1/6];
        
        % run the simulation
        options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3,'Events',detect_blowup);
        [t3,u_raw3] = ode45(@(t,u) RHS(u,t,params3),[0,end_time],u_init(:),options);
        
        
        % reshape the output array into an intelligible shape (should make this a
        % separate function later)
        u_array3 = zeros([size(u_init) length(t3)]);
        for i = 1:length(t3)
            u_array3(:,:,:,:,i) = reshape(u_raw3(i,:),[N,N,2,2]);
        end
        
        save(sprintf('t3_%i_%i',N,end_time),'t3');
        save(sprintf('u_array3_%i_%i.mat',N,end_time),'u_array3');
    
    params2 = params;
    params2.func = @(x) t2model_RHS(x);
    params2.coeff = [1 -1/2];
        
        % run the simulation
        options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3,'Events',detect_blowup);
        [t2,u_raw2] = ode45(@(t,u) RHS(u,t,params2),[0,end_time],u_init(:),options);
        
        
        % reshape the output array into an intelligible shape (should make this a
        % separate function later)
        u_array2 = zeros([size(u_init) length(t2)]);
        for i = 1:length(t2)
            u_array2(:,:,:,:,i) = reshape(u_raw2(i,:),[N,N,2,2]);
        end
        
        save(sprintf('t2_%i_%i',N,end_time),'t2');
        save(sprintf('u_array2_%i_%i.mat',N,end_time),'u_array2');
        
    params1 = params;
    params1.func = @(x) tmodel_RHS(x);
    params1.coeff = 1;
        
        % run the simulation
        options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3,'Events',detect_blowup);
        [t1,u_raw1] = ode45(@(t,u) RHS(u,t,params1),[0,end_time],u_init(:),options);
        
        
        % reshape the output array into an intelligible shape (should make this a
        % separate function later)
        u_array1 = zeros([size(u_init) length(t1)]);
        for i = 1:length(t1)
            u_array1(:,:,:,:,i) = reshape(u_raw1(i,:),[N,N,2,2]);
        end
        
        save(sprintf('t1_%i_%i',N,end_time),'t1');
        save(sprintf('u_array1_%i_%i.mat',N,end_time),'u_array1');
  
        
        energy_exact = get_2D_energy(u,N);
        energy_t = get_2D_energy(u_array1,N);
        energy_t2 = get_2D_energy(u_array2,N);
        energy_t3 = get_2D_energy(u_array3,N);
        energy_t4 = get_2D_energy(u_array4,N);
        
        figure
        hold on
        plot(log(t),log(energy_exact))
        plot(log(t1),log(energy_t))
        plot(log(t2),log(energy_t2))
        plot(log(t3),log(energy_t3))
        plot(log(t4),log(energy_t4))
        
        title(sprintf('Energy in resolved modes, N=%i',N))
        legend('exact','t-model','t^2-model','t^3-model','t^4-model');
        xlabel('log(time)')
        ylabel('log(energy)')
end
    
