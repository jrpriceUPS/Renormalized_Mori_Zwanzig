clear all;close all;

format long

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

end_time = 3.3738; % point of instability

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

params_t3model.k = k;
params_t3model.N = N;
params_t3model.M = M;
params_t3model.func = @(x) t3model_RHS(x);
params_t3model.coeff = [1,1,1];
params_t3model.a = 2:M;
params_t3model.b = 2*M:-1:M+2;
params_t3model.a_tilde = N+1:M;
params_t3model.print_time = 1;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on');
[t_t3model,u_raw_t3model] = ode45(@(t,u) RHS(u,t,params_t3model),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_t3model = zeros([size(u) length(t_t3model)]);
for i = 1:length(t_t3model)
    u_array_t3model(:,:,:,:,:,i) = reshape(u_raw_t3model(i,:),[N,N,N,3,4]);
end

save t_t3model t_t3model
save u_array_t3model u_array_t3model

load('t_full');
load('u_array_full');

load('t_markov')
load('u_array_markov')

load('t_tmodel')
load('u_array_tmodel')

load('t_t2model')
load('u_array_t2model')

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N);
energy_tmodel = get_3D_energy(u_array_tmodel,N);
energy_markov = get_3D_energy(u_array_markov,N);
energy_t2model = get_3D_energy(u_array_t2model,N);
energy_t3model = get_3D_energy(u_array_t3model,N);
figure(1)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_markov,energy_markov,'r','linewidth',2)
plot(t_tmodel,energy_tmodel,'k','linewidth',2)
plot(t_t2model,energy_t2model,'g','linewidth',2)
plot(t_t3model,energy_t3model,'c','linewidth',2)
axis([0,8,0,0.5])
legend('Full model, N = 8','Markov model, N = 4','t-model, N = 4','t^2-model, N = 4','t^3-model, N = 4','location','southwest')
title('Energy in first N = 4 modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy1','png')

% plot the energy in some modes
energy_full = get_3D_energy(u_array_full,N/2);
energy_tmodel = get_3D_energy(u_array_tmodel,N/2);
energy_markov = get_3D_energy(u_array_markov,N/2);
energy_t2model = get_3D_energy(u_array_t2model,N/2);
energy_t3model = get_3D_energy(u_array_t3model,N/2);
figure(2)
hold off
plot(t_full,energy_full,'linewidth',2)
hold on
plot(t_markov,energy_markov,'r','linewidth',2)
plot(t_tmodel,energy_tmodel,'k','linewidth',2)
plot(t_t2model,energy_t2model,'g','linewidth',2)
plot(t_t3model,energy_t3model,'c','linewidth',2)
axis([0,8,0,0.5])
legend('Full model, N = 8','Markov model, N = 4','t-model, N = 4','t^2-model, N = 4','t^3-model, N = 4','location','southwest')
title('Energy in first N = 2 modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy2','png')

% plot the helicity
w_full = helicity(u_array_full);
w_tmodel = helicity(u_array_tmodel);
w_markov = helicity(u_array_markov);
w_t2model = helicity(u_array_t2model);
w_t3model = helicity(u_array_t3model);
figure(3)
plot(t_full,w_full,'linewidth',2)
hold on
plot(t_markov,w_markov,'r','linewidth',2)
plot(t_tmodel,w_tmodel,'k','linewidth',2)
plot(t_t2model,w_t2model,'g','linewidth',2)
plot(t_t3model,w_t3model,'c','linewidth',2)
legend('Full model, N = 8','Markov model, N = 4','t-model, N = 4','t^2-model, N = 4','t^3-model, N = 4','location','southwest')
title('Helicity','fontsize',16)
xlabel('time','fontsize',16)
ylabel('w','fontsize',16)
saveas(gcf,'helicity','png')