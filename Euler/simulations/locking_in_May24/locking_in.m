clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N = 48;
M = 3*N;


if ~(exist('u48.mat','file') == 2)
    create_data(N,1);
end
if ~(exist('u48_2.mat','file')==2)
    
    load u48;
    load(sprintf('u%i.mat',N))
    load(sprintf('t%i.mat',N))
    start_time = t(end);
    u0 = u(:,:,:,:,:,end);
    
    
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
    options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
    [t_new,u_raw] = ode45(@(t,u) RHS(u,t,params),[1,1.5],u0(:),options);
    
    % reshape the output array into an intelligible shape (should make this a
    % separate function later)
    u_new = zeros([size(u0) length(t_new)]);
    for i = 1:length(t_new)
        u_new(:,:,:,:,:,i) = reshape(u_raw(i,:),[N,N,N,3,4]);
    end
    
    t2 = t_new;
    u2 = u_new;
    
    save(sprintf('t%i_2',N),'t2');
    save(sprintf('u%i_2',N),'u2');
    
end




load u48
load u48_2

load t48
load t48_2

s = size(u);

u_both = zeros(s(1),s(2),s(3),s(4),s(5),length(t)+length(t2));
u_both(:,:,:,:,:,1:length(t)) = u;
u_both(:,:,:,:,:,length(t)+1:end) = u2;

t_both = [t;t2];

if ~(exist('tmodel_size_list48.mat','file') == 2)
    [tmodel_size_list,tmodel_size_list_full] = resolve_array(u_both,t_both);
    save('tmodel_size_list48','tmodel_size_list')
    save('tmodel_size_list_full48','tmodel_size_list_full')
end

load tmodel_size_list48
min_tol = 1e-16;
max_tol = 1e-10;

viable_snapshots = find(tmodel_size_list > min_tol & tmodel_size_list < max_tol);
snaps = 10:5:length(viable_snapshots);

N_list = 4:2:22;

coeffs = zeros(4,length(N_list),4,length(snaps));
coeffs_t = zeros(4,length(N_list),4,length(snaps));


for i = 1:length(snaps)
    window = viable_snapshots(1:snaps(i));
    
    % trim the arrays to those viable times
    u_array = u_both(:,:,:,:,:,window);
    t_array = t_both(window);
    
    % compute the renormalization coefficients
    [c,~,~] = renormalize(u_array,N_list,t_array,0,1);
    coeffs(:,:,:,i) = c;
    
    [c_t,~,~] = renormalize(u_array,N_list,t_array,1,1);
    coeffs_t(:,:,:,i) = c_t;
    
end

figure(1)
hold off
plot(t_array(snaps),squeeze(coeffs(1,:,1,:)))
title('Only fit R1 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R1','png')




figure(2)
hold off
subplot(2,1,1)
plot(t_array(snaps),squeeze(coeffs(1,:,2,:)))
title('R1+R2 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,1,2)
plot(t_array(snaps),squeeze(coeffs(2,:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R2','png')



figure(3)
hold off
subplot(3,1,1)
plot(t_array(snaps),squeeze(coeffs(1,:,3,:)))
title('R1+R2+R3 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,2)
plot(t_array(snaps),squeeze(coeffs(2,:,3,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,3)
plot(t_array(snaps),squeeze(coeffs(3,:,3,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R3','png')



figure(4)
hold off
subplot(2,2,1)
plot(t_array(snaps),squeeze(coeffs(1,:,4,:)))
title('R1+R2+R3+R4 model, decaying','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,2)
plot(t_array(snaps),squeeze(coeffs(2,:,4,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,3)
plot(t_array(snaps),squeeze(coeffs(3,:,4,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,4)
plot(t_array(snaps),squeeze(coeffs(4,:,4,:)))
ylabel('R4','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R4','png')










figure(5)
hold off
plot(t_array(snaps),squeeze(coeffs_t(1,:,1,:)))
title('Only fit R1 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R1_c','png')




figure(6)
hold off
subplot(2,1,1)
plot(t_array(snaps),squeeze(coeffs_t(1,:,2,:)))
title('R1+R2 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,1,2)
plot(t_array(snaps),squeeze(coeffs_t(2,:,2,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R2_c','png')



figure(7)
hold off
subplot(3,1,1)
plot(t_array(snaps),squeeze(coeffs_t(1,:,3,:)))
title('R1+R2+R3 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,2)
plot(t_array(snaps),squeeze(coeffs_t(2,:,3,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(3,1,3)
plot(t_array(snaps),squeeze(coeffs_t(3,:,3,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R3_c','png')



figure(8)
hold off
subplot(2,2,1)
plot(t_array(snaps),squeeze(coeffs_t(1,:,4,:)))
title('R1+R2+R3+R4 model, constant','fontsize',16)
ylabel('R1','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,2)
plot(t_array(snaps),squeeze(coeffs_t(2,:,4,:)))
ylabel('R2','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,3)
plot(t_array(snaps),squeeze(coeffs_t(3,:,4,:)))
ylabel('R3','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
subplot(2,2,4)
plot(t_array(snaps),squeeze(coeffs_t(4,:,4,:)))
ylabel('R4','fontsize',16)
xlabel('Edge of window T ([0,T])','fontsize',16)
saveas(gcf,'R4_c','png')
