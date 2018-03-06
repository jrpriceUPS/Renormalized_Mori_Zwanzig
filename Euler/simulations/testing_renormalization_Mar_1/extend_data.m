

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

load t_save
load u_save

N = 16;
M = 3*N;

start_time = t(end);
initial_condition = u_save(:,:,:,:,:,end);

u0 = initial_condition(:);

end_time = start_time + 1;

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
[t_full,u_raw_full] = ode45(@(t,u) RHS(u,t,params),[start_time,end_time],u0(:),options);

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_full_data = zeros([size(u) length(t_full)]);
for i = 1:length(t_full)
    u_full_data(:,:,:,:,:,i) = reshape(u_raw_full(i,:),[N,N,N,3,4]);
end

t = [t;t_full(2:end)];

s = size(u);

new_u = zeros(s(1),s(2),s(3),s(4),s(5),length(t));
new_u(:,:,:,:,:,1:start_time) = u;
new_u(:,:,:,:,:,start_time+1:end) = u_full_data(:,:,:,:,:,2:end);
u = new_u;

save u u
save t t
