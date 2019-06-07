function create_data_alpha(N,end_time,alpha,extension)
%
% A function to generate a full model simulation of size N up to time
% end_time
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%         N  =  resolution of full model  
%
%  end_time  =  final time desired for simulation
%
%     alpha  =  size of the Taylor-Green initial condition
%
% extension  =  extension to be appended to name of output to differentiate
%
%
%%%%%%%%%
%OUTPUS:%
%%%%%%%%%
%
%  uN.dat  =  saved array of Fourier modes from time zero to time end_time
%             (N x N x N x 3 x 4 x length(tN))
%
%  tN.dat  =  saved array of times associated with solution

% load relevant folders into the path
addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

% size of array needed for dealiasing
M = 3*N;

if exist(['u48_' extension '.mat'],'file') == 2
    
    % if there is already data for this resolution, load it and continue the
    % simulation from the previous end time up to the proposed end time
    load(['u48_' extension '.mat'])
    load(['t48_' extension '.mat'])
    start_time = t(end);
    u0 = u(:,:,:,:,:,end);
    
else
    
    % otherwise, start at time 0 and use Taylor-Green for the initial
    % condition
    start_time = 0;
    
    % uniform grid
    x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
    y = x;
    z = x;
    
    % 3D array of data points
    [X,Y,Z] = ndgrid(x,y,z);
    
    % create initial condition
    eval = alpha*taylor_green(X,Y,Z);
    u_full = fftn_norm(eval);
    u0 = u_squishify(u_full,N);
    
end

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

% load relevant parameters into parameter structure
params.k = k;
params.N = N;
params.M = M;
params.func = @(x) full_RHS(x);
params.coeff = [];
params.a = 2:M;
params.b = 2*M:-1:M+2;
params.a_tilde = N+1:M;
params.a_tilde2 = 2*N+1:M;
params.print_time = 1;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t_new,u_raw] = ode45(@(t,u) RHS(u,t,params),[start_time,end_time],u0(:),options);

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_new = zeros([size(u0) length(t_new)]);
for i = 1:length(t_new)
    u_new(:,:,:,:,:,i) = reshape(u_raw(i,:),[N,N,N,3,4]);
end

if exist(['u48_' extension '.mat'],'file') == 2
    
    % if this resolution has been run before, append the new results to
    % those results
    u(:,:,:,:,:,length(t)+1:length(t)+length(t_new)-1) = u_new(:,:,:,:,:,2:end);
    t = [t;t_new(2:end)];
    
else
    
    % otherwise, save the current results
    u = u_new;
    t = t_new;
    
end

% save the results into the directory
save(['u48_' extension '.mat'],'t');
save(['t48_' extension '.mat'],'u');