clear all;close all;

format long

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

end_time = 100; % point of instability

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

params_t4model.k = k;
params_t4model.N = N;
params_t4model.M = M;
params_t4model.func = @(x) t4model_RHS(x);
params_t4model.coeff = [0.584216814927850,0.278260135212772,0.084385844450150,0.010981993693650];
params_t4model.a = 2:M;
params_t4model.b = 2*M:-1:M+2;
params_t4model.a_tilde = N+1:M;
params_t4model.print_time = 1;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t_t4model,u_raw_t4model] = ode45(@(t,u) RHS(u,t,params_t4model),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array_t4model = zeros([size(u) length(t_t4model)]);
for i = 1:length(t_t4model)
    u_array_t4model(:,:,:,:,:,i) = reshape(u_raw_t4model(i,:),[N,N,N,3,4]);
end

save t_t4model t_t4model
save u_array_t4model u_array_t4model

% plot the energy in some modes
energy_t4model = get_3D_energy(u_array_t4model,N);
figure(1)
hold off
plot(log(t_t4model),log(energy_t4model),'linewidth',2)
legend('Renormalized model, N = 4','location','southwest')
title('Energy in resolved modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'energy','png')



% plot the helicity
w_t4model = helicity(u_array_t4model);
figure(2)
hold off
plot(t_t4model,w_t4model,'m','linewidth',2)
legend('Renormalized model, N = 4','location','southwest')
title('Helicity','fontsize',16)
xlabel('time','fontsize',16)
ylabel('w','fontsize',16)
saveas(gcf,'helicity','png')