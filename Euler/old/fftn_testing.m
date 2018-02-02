% A simulation script to run a full 3D Euler simulation and make sure
% convolution functions and analysis functions are functioning properly


clear all;close all;

addpath ././nonlinear
addpath ././simulation_functions
addpath ././analysis

N = 4; % resolution
M = 3*N; % size needed for dealising full model
end_time = 1; % end time of simulation

% uniform grid
x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

% taylor-green initial condition
eval = taylor_green(X,Y,Z);

% FFT to get initial condition
u_full = fftn_norm(eval);

% reshape into compact form
u = u_squishify(u_full,N);

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on');
[t,u_raw] = ode45(@(t,u) full_euler(u,k,N,M),[0,end_time],u(:),options);

% reshape the array into the most convenient shape
u_array = zeros([size(u) length(t)]);
for i = 1:length(t)
    u_array(:,:,:,:,:,i) = reshape(u_raw(i,:),[N,N,N,3,4]);
end

% plot total energy (should be constant)
energy = get_3D_energy(u_array,N);
figure(1)
hold off
plot(t,energy,'linewidth',2)
title('Total energy','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'total_energy','png')

% plot energy in a subset of the modes (should not be monotonic)
energy = get_3D_energy(u_array,N/2);
figure(2)
hold off
plot(t,energy,'linewidth',2)
title('Energy in half the modes','fontsize',16)
xlabel('time','fontsize',16)
ylabel('energy','fontsize',16)
saveas(gcf,'subset_energy','png')

% plot the helicity (should be zero)
w = helicity(u_array);
figure(3)
plot(t,w,'linewidth',2)
title('Helicity','fontsize',16)
xlabel('time','fontsize',16)
ylabel('w','fontsize',16)
saveas(gcf,'helicity','png')