clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

N_full = 48;
alpha = 1;
N_list_coeff = 4:2:24;
tau_list = [0.0;0.2;0.4;0.6;0.8;1.0];
N_list = 6:2:14;

% Load calculated coefficients
% This data is generated using c_vs_N_test_coeff_data.m
load(['c1_op_n_1_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c2_op_n_2_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c3_op_n_3_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list_coeff(1)) '_to_' num2str(N_list_coeff(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])

% Choose tau value
tau = 1.0;

% Find coefficients corresponding to tau and N
ind = find(tau == tau_list);
c1_op_tau = c1_op{1,ind};
c2_op_tau = c2_op{1,ind};
c3_op_tau = c3_op{1,ind};
c4_op_tau = c4_op{1,ind};

ind_N_l = find(N_list_coeff == N_list(1));
ind_N_h = find(N_list_coeff == N_list(end));
c1_op_N = c1_op_tau(:,ind_N_l:ind_N_h);
c2_op_N = c2_op_tau(:,ind_N_l:ind_N_h);
c3_op_N = c3_op_tau(:,ind_N_l:ind_N_h);
c4_op_N = c4_op_tau(:,ind_N_l:ind_N_h);

[c1_laws,r1] = create_scaling_laws(N_list(1:end),c1_op_N(1,:));
[c2_laws,r2] = create_scaling_laws(N_list(1:end),c2_op_N(1:2,:));
[c3_laws,r3] = create_scaling_laws(N_list(1:end),c3_op_N(1:3,:));
[c4_laws,r4] = create_scaling_laws(N_list(1:end),c4_op_N(1:4,:));

% Set N for simulation
N = 16;
sim_endtime = 1000;

M = 3*N;

% Calculate coefficients from scaling laws
coeffs = zeros(5,1);
    
coeffs(1) = exp(c4_laws(1,2))*N^c4_laws(1,1);
coeffs(2) = -exp(c4_laws(2,2))*N^c4_laws(2,1);
coeffs(3) = exp(c4_laws(3,2))*N^c4_laws(3,1);
coeffs(4) = -exp(c4_laws(4,2))*N^c4_laws(4,1);
coeffs(5) = tau;

% Uniform grid
x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

% Create initial condition
eval = taylor_green(X,Y,Z)*alpha;
u_full = fftn_norm(eval);
u = u_squishify(u_full,N);

% k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

params.k = k;
params.N = N;
params.M = M;
params.a = 2:M;
params.b = 2*M:-1:M+2;
params.a_tilde = N+1:M;
params.a_tilde2 = 2*N+1:M;
params.print_time = 1;
params.tau = coeffs(5);
params4 = params;
params4.func = @(x) t4model_RHS_tau(x);
params4.coeff = coeffs(1:4);

% Isolate modes that do not include any components larger than N
a = [1:N,2*M-N+1:2*M];

percent_increase = 10;
u_full = u_fullify(u,M);
init_energy = sum(sum(sum(sum(u_full(a,a,a,:).*conj(u_full(a,a,a,:))))));
options = odeset('AbsTol',1e-14,'Stats','on','InitialStep',1e-3,'Events',@(t,u) energy_grow(t,u,init_energy,percent_increase,params4));
[t,u_raw] = ode45(@(t,u) RHS(u,t,params4),[0,sim_endtime],u(:),options);

% Reshape the output array into an intelligible shape
u = zeros([size(u) length(t)]);
for j = 1:length(t)
    u(:,:,:,:,:,j) = reshape(u_raw(j,:),[N,N,N,3,4]);
end

% Calculate the energy in the N modes of the ROM
energy = get_3D_energy(u,N);

% Save the results into the directory
save(['u_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'u','-v7.3');
save(['t_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'t');
