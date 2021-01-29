function [u,t] = create_data(N,end_time,alpha)
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
%%%%%%%%%
%OUTPUS:%
%%%%%%%%%
%
%  u_N_endtime_x10_end_time.mat  =  saved array of Fourier modes from time zero to time end_time
%             (N x N x N x 3 x 4 x length(tN))
%
%  t_N_endtime_x10_end_time.mat  =  saved array of times associated with solution

% load relevant folders into the path
addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

% size of array needed for dealiasing
M = 3*N;

%Taylor-Green for the initial condition
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
options = odeset('AbsTol',1e-14,'Stats','on','InitialStep',1e-3);
[t,u_raw] = ode45(@(t,u) RHS(u,t,params),[start_time,end_time],u0(:),options);

% reshape the output array into an intelligible shape
u = zeros([size(u0) length(t)]);
for i = 1:length(t)
    u(:,:,:,:,:,i) = reshape(u_raw(i,:),[N,N,N,3,4]);
end

% save the results into the directory
save(['u_' num2str(N) '_endtime_x10_' num2str(end_time*10) '.mat'],'u','-v7.3');
save(['t_' num2str(N) '_endtime_x10_' num2str(end_time*10) '.mat'],'t');