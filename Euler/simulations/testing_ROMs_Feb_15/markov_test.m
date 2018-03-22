% original testbed for full Euler implementation

clear all;close all;
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

end_time = 8;

N_full = 8; % resolution
M_full = 3*N_full;

N_markov = 4;
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
params_full.func = @(x) full_RHS(x);
params_full.coeff = [];
params_full.a = 2:M_full;
params_full.b = 2*M_full:-1:M_full+2;
params_full.a_tilde = N_full+1:M_full;
params_full.print_time = 1;

params_markov.k = k_markov;
params_markov.N = N_markov;
params_markov.M = M_markov;
params_markov.func = @(x) markov_RHS(x);
params_markov.coeff = [];
params_markov.a = 2:M_markov;
params_markov.b = 2*M_markov:-1:M_markov+2;
params_markov.a_tilde = N_markov+1:M_markov;
params_markov.a_tilde2 = 2*N_markov+1:M_markov;
params_markov.print_time = 1;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t_full,u_raw_full] = ode45(@(t,u) RHS(u,t,params_full),[0,end_time],u1(:),options);
[t_markov,u_raw_markov] = ode45(@(t,u) RHS(u,t,params_markov),[0,end_time],u2(:),options);

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_full = zeros([size(u1) length(t_full)]);
for i = 1:length(t_full)
    u_array_full(:,:,:,:,:,i) = reshape(u_raw_full(i,:),[N_full,N_full,N_full,3,4]);
end

% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_markov = zeros([size(u2) length(t_markov)]);
for i = 1:length(t_markov)
    u_array_markov(:,:,:,:,:,i) = reshape(u_raw_markov(i,:),[N_markov,N_markov,N_markov,3,4]);
end

save u_array_full u_array_full
save u_array_markov u_array_markov

save t_full t_full
save t_markov t_markov

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N_markov);
energy_markov = get_3D_energy(u_array_markov,N_markov);
figure(1)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_markov,energy_markov,'r','linewidth',2)
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