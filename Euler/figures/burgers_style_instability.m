function burgers_style_instability(N,end_time,filetype)

format long

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

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

params.k = k;
params.N = N;
params.M = M;
params.a = 2:M;
params.b = 2*M:-1:M+2;
params.a_tilde = N+1:M;
params.a_tilde2 = 2*N+1:M;
params.print_time = 1;
params.no_time = 0;

params0 = params;
params0.func = @(x) markov_RHS(x);

if exist(sprintf('u_array0_%i_%i.mat',N,end_time),'file') == 2
    
    load(sprintf('u_array0_%i_%i.mat',N,end_time))
    load(sprintf('t0_%i_%i',N,end_time))
    
else
    
    
    % run the simulation
    options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
    [t0,u_raw0] = ode45(@(t,u) RHS(u,t,params0),[0,end_time],u(:),options);
    
    
    
    
    % reshape the output array into an intelligible shape (should make this a
    % separate function later)
    u_array0 = zeros([size(u) length(t0)]);
    for j = 1:length(t0)
        u_array0(:,:,:,:,:,j) = reshape(u_raw0(j,:),[N,N,N,3,4]);
    end
    
    save(sprintf('t0_%i_%i',N,end_time),'t0');
    save(sprintf('u_array0_%i_%i',N,end_time),'u_array0');
    
end

params1 = params;
params1.func = @(x) tmodel_RHS(x);
params1.coeff = scaling_law_time(N,1);


% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t1,u_raw1] = ode45(@(t,u) RHS(u,t,params1),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array1 = zeros([size(u) length(t1)]);
for i = 1:length(t1)
    u_array1(:,:,:,:,:,i) = reshape(u_raw1(i,:),[N,N,N,3,4]);
end

save(sprintf('t1_%i_%it',N,end_time),'t1');
save(sprintf('u_array1_%i_%it',N,end_time),'u_array1');


params2 = params;
params2.func = @(x) t2model_RHS(x);
params2.coeff = scaling_law(N,2);

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t2,u_raw2] = ode45(@(t,u) RHS(u,t,params2),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array2 = zeros([size(u) length(t2)]);
for i = 1:length(t2)
    u_array2(:,:,:,:,:,i) = reshape(u_raw2(i,:),[N,N,N,3,4]);
end

save(sprintf('t2_%i_%it',N,end_time),'t2');
save(sprintf('u_array2_%i_%it',N,end_time),'u_array2');



params3 = params;
params3.func = @(x) t3model_RHS(x);
params3.coeff = scaling_law(N,3);

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t3,u_raw3] = ode45(@(t,u) RHS(u,t,params3),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array3 = zeros([size(u) length(t3)]);
for i = 1:length(t3)
    u_array3(:,:,:,:,:,i) = reshape(u_raw3(i,:),[N,N,N,3,4]);
end

save(sprintf('t3_%i_%it',N,end_time),'t3');
save(sprintf('u_array3_%i_%it',N,end_time),'u_array3');

params4 = params;
params4.func = @(x) t4model_RHS(x);
params4.coeff = scaling_law(N,4);

% run the simulation
options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
[t4,u_raw4] = ode45(@(t,u) RHS(u,t,params4),[0,end_time],u(:),options);


% reshape the output array into an intelligible shape (should make this a
% separate function later)
u_array4 = zeros([size(u) length(t4)]);
for i = 1:length(t4)
    u_array4(:,:,:,:,:,i) = reshape(u_raw4(i,:),[N,N,N,3,4]);
end

save(sprintf('t4_%i_%it',N,end_time),'t4');
save(sprintf('u_array4_%i_%it',N,end_time),'u_array4');




% plot the energy in some modes
energy0 = get_3D_energy(u_array0,N);
energy1 = get_3D_energy(u_array1,N);
energy2 = get_3D_energy(u_array2,N);
energy3 = get_3D_energy(u_array3,N);
energy4 = get_3D_energy(u_array4,N);
figure(1)
hold off
plot(log(t0),log(energy0),'b','linewidth',2)
hold on
plot(log(t1),log(energy1),'r','linewidth',2)
plot(log(t2),log(energy2),'g','linewidth',2)
plot(log(t3),log(energy3),'k','linewidth',2)
plot(log(t4),log(energy4),'c','linewidth',2)
legend('Markov ROM','1st order ROM','2nd order ROM','3rd order ROM','4th order ROM','location','southwest')
title(sprintf('Energy in resolved modes, N = %i',N),'fontsize',16)
xlabel('log(time)','fontsize',16)
ylabel('log(energy)','fontsize',16)
saveas(gcf,sprintf('unstable_burgers_%i',N),filetype)




