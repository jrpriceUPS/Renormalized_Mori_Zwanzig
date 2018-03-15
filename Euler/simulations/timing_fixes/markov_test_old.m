% original testbed for full Euler implementation

clear all;close all;
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis
addpath ../timing_fixes

end_time = 8;

N_full = 4; % resolution
M_full = 3*N_full;

N_markov = 2;
M_markov = 3*N_markov;

% uniform grid
x = linspace(0,2*pi*(2*M_full-1)/(2*M_full),2*M_full).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

% create initial condition
eval = taylor_green(X,Y,Z);
u_full = fftn_norm(eval);
u1 = u_squishify(u_full,N_full);
u2 = u_squishify(u_full,N_markov);

% make k array
k_vec_full = [0:M_full-1,-M_full:1:-1];
[kx_full,ky_full,kz_full] = ndgrid(k_vec_full,k_vec_full,k_vec_full);
k_full = zeros(2*M_full,2*M_full,2*M_full,3);
k_full(:,:,:,1) = kx_full;
k_full(:,:,:,2) = ky_full;
k_full(:,:,:,3) = kz_full;

% make k array
k_vec_markov = [0:M_markov-1,-M_markov:1:-1];
[kx_markov,ky_markov,kz_markov] = ndgrid(k_vec_markov,k_vec_markov,k_vec_markov);
k_markov = zeros(2*M_markov,2*M_markov,2*M_markov,3);
k_markov(:,:,:,1) = kx_markov;
k_markov(:,:,:,2) = ky_markov;
k_markov(:,:,:,3) = kz_markov;

params_full.k = k_full;
params_full.N = N_full;
params_full.M = M_full;
params_full.func = @(x) full_RHS_old(x);
params_full.coeff = [];
params_full.a = 2:M_full;
params_full.b = 2*M_full:-1:M_full+2;
params_full.a_tilde = N_full+1:M_full;
params_full.print_time = 1;

params_markov.k = k_markov;
params_markov.N = N_markov;
params_markov.M = M_markov;
params_markov.func = @(x) markov_RHS_old(x);
params_markov.coeff = [];
params_markov.a = 2:M_markov;
params_markov.b = 2*M_markov:-1:M_markov+2;
params_markov.a_tilde = N_markov+1:M_markov;
params_markov.print_time = 1;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t_full_old,u_raw_full_old] = ode45(@(t,u) RHS(u,t,params_full),[0,end_time],u1(:),options);
[t_markov_old,u_raw_markov_old] = ode45(@(t,u) RHS(u,t,params_markov),[0,end_time],u2(:),options);

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_full_old = zeros([size(u1) length(t_full_old)]);
for i = 1:length(t_full_old)
    u_array_full_old(:,:,:,:,:,i) = reshape(u_raw_full_old(i,:),[N_full,N_full,N_full,3,4]);
end

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_markov_old = zeros([size(u2) length(t_markov_old)]);
for i = 1:length(t_markov_old)
    u_array_markov_old(:,:,:,:,:,i) = reshape(u_raw_markov_old(i,:),[N_markov,N_markov,N_markov,3,4]);
end

save u_array_full_old u_array_full_old
save u_array_markov_old u_array_markov_old

save t_full_old t_full_old
save t_markov_old t_markov_old

load u_array_full
load u_array_markov
load t_full
load t_markov

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N_markov);
energy_markov = get_3D_energy(u_array_markov,N_markov);
energy_full_old = get_3D_energy(u_array_full_old,N_markov);
energy_markov_old = get_3D_energy(u_array_markov_old,N_markov);
figure(1)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_markov,energy_markov,'r','linewidth',2)
plot(t_full_old,energy_full_old,'g','linewidth',2)
plot(t_markov_old,energy_markov_old,'k','linewidth',2)
legend('Full model, N = 8','Markov model, N = 4')
title('Energy in first N = 4 modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy1','png')

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N_markov/2);
energy_markov = get_3D_energy(u_array_markov,N_markov/2);
figure(2)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_markov,energy_markov,'r','linewidth',2)
legend('Full model, N = 8','Markov model, N = 4')
title('Energy in first N = 2 modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy2','png')

% plot the helicity
w_full = helicity(u_array_full);
w_markov = helicity(u_array_markov);
figure(3)
plot(t_full,w_full,'linewidth',2)
hold on
plot(t_markov,w_markov,'r','linewidth',2)
title('Helicity','fontsize',16)
xlabel('time','fontsize',16)
ylabel('w','fontsize',16)
saveas(gcf,'helicity','png')