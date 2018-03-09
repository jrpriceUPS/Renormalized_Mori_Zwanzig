function create_data(N,end_time)

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

M = 3*N;

if exist(sprintf('u%i.mat',N),'file') == 2
    
    load(sprintf('u%i.mat',N))
    load(sprintf('t%i.mat',N))
    start_time = t(end);
    u0 = u(:,:,:,:,:,end);
    
else
    
    start_time = 0;
    
    % uniform grid
    x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
    y = x;
    z = x;
    
    % 3D array of data points
    [X,Y,Z] = ndgrid(x,y,z);
    
    % create initial condition
    eval = taylor_green(X,Y,Z);
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


params.k = k;
params.N = N;
params.M = M;
params.func = @(x) full_RHS(x);
params.coeff = [];
params.a = 2:M;
params.b = 2*M:-1:M+2;
params.a_tilde = N+1:M;
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

if exist(sprintf('u%i.mat',N),'file') == 2
    
    u(:,:,:,:,:,length(t)+1:length(t)+length(t_new)-1) = u_new(:,:,:,:,:,2:end);
    t = [t;t_new(2:end)];
    
    
    
else
    
    u = u_new;
    t = t_new;
    
end

save(sprintf('t%i',N),'t');
save(sprintf('u%i',N),'u');