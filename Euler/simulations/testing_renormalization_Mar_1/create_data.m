clear all;close all;
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

end_time = 2;

N = 16; % resolution
M = 3*N;

% uniform grid
x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

% create initial condition
eval = taylor_green(X,Y,Z);
u_full = fftn_norm(eval);
u = u_squishify(u_full,N);

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
options = odeset('RelTol',1e-10,'Stats','on');
[t_full,u_raw_full] = ode45(@(t,u) RHS(u,t,params),[0,end_time],u(:),options);

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_full_data = zeros([size(u) length(t_full)]);
for i = 1:length(t_full)
    u_full_data(:,:,:,:,:,i) = reshape(u_raw_full(i,:),[N,N,N,3,4]);
end

save t_full t_full

save u_full_data u_full_data