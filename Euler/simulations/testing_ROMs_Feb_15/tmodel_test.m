%clear all;close all;

format long

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

end_time = 8;

N = 4;
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

params_tmodel.k = k;
params_tmodel.N = N;
params_tmodel.M = M;
params_tmodel.func = @(x) tmodel_RHS(x);
params_tmodel.coeff = 1;
params_tmodel.a = 2:M;
params_tmodel.b = 2*M:-1:M+2;
params_tmodel.a_tilde = N+1:M;
params_tmodel.a_tilde2 = 2*N+1:M;
params_tmodel.print_time = 1;
params_tmodel.no_time = 0;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on');
[t_tmodel,u_raw_tmodel] = ode45(@(t,u) RHS(u,t,params_tmodel),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_tmodel = zeros([size(u) length(t_tmodel)]);
for i = 1:length(t_tmodel)
    u_array_tmodel(:,:,:,:,:,i) = reshape(u_raw_tmodel(i,:),[N,N,N,3,4]);
end

save t_tmodel t_tmodel
save u_array_tmodel u_array_tmodel

load('t_full');
load('u_array_full');

load('t_markov')
load('u_array_markov')

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N);
energy_tmodel = get_3D_energy(u_array_tmodel,N);
energy_markov = get_3D_energy(u_array_markov,N);
figure(1)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_tmodel,energy_tmodel,'r','linewidth',2)
plot(t_markov,energy_markov,'k','linewidth',2)
legend('Full model, N = 8','t-model, N = 4','Markov model, N = 4')
title('Energy in first N = 4 modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy1','png')

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N/2);
energy_tmodel = get_3D_energy(u_array_tmodel,N/2);
energy_markov = get_3D_energy(u_array_markov,N/2);
figure(2)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_tmodel,energy_tmodel,'r','linewidth',2)
plot(t_markov,energy_markov,'k','linewidth',2)
legend('Full model, N = 8','t-model, N = 4','Markov model, N = 4')
title('Energy in first N = 2 modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy2','png')

% plot the helicity
w_full = helicity(u_array_full);
w_tmodel = helicity(u_array_tmodel);
w_markov = helicity(u_array_markov);
figure(3)
plot(t_full,w_full,'linewidth',2)
hold on
plot(t_tmodel,w_tmodel,'r','linewidth',2)
plot(t_markov,w_markov,'k','linewidth',2)
legend('Full model, N = 8','t-model, N = 4','Markov model, N = 4')
title('Helicity','fontsize',16)
xlabel('time','fontsize',16)
ylabel('w','fontsize',16)
saveas(gcf,'helicity','png')