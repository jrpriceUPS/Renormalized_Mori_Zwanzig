clear all; close all; clc;

addpath ../../simulation_functions/
addpath ../../nonlinear/
addpath ../../analysis/

N_full = 48;
sim_endtime = 1000;
alpha = 1;
N_list = 4:2:24;
N = 6;
tau_list = [0.0;0.2;0.4;0.6;0.8;1.0];
tau = 1.0;
c_list = (4/5)*hsv(length(tau_list));

% Load calculated coefficients
% This data is generated using c_vs_N_test_coeff_data.m
load(['c1_op_n_1_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c2_op_n_2_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c3_op_n_3_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])
load(['c4_op_n_4_M_' num2str(N_full) '_N_' num2str(N_list(1)) '_to_' num2str(N_list(end)) '_tau_' num2str(tau_list(1)) '_to_' num2str(tau_list(end)) '.mat'])

% t-model

figure()
hold on
for i = 1:length(tau_list)
    
    tau = tau_list(i);
    
    % Find coefficients corresponding to tau and N
    ind = find(tau == tau_list); 
    c1_op_tau = c1_op{1,ind};
    ind_N = find(N == N_list); 
    c1_op_N_tau = c1_op_tau(:,ind_N);

    M = 3*N;
    
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
    params.tau = c1_op_N_tau(2);
    params1 = params;
    params1.func = @(x) tmodel_RHS_tau(x);
    params1.coeff = c1_op_N_tau(1);
    
    % Isolate modes that do not include any components larger than N
    a = [1:N,2*M-N+1:2*M];
    
    percent_increase = 10;
    u_full = u_fullify(u,M);
    init_energy = sum(sum(sum(sum(u_full(a,a,a,:).*conj(u_full(a,a,a,:))))));
    options = odeset('AbsTol',1e-14,'Stats','on','InitialStep',1e-3,'Events',@(t,u) energy_grow(t,u,init_energy,percent_increase,params1));
    [t,u_raw] = ode45(@(t,u) RHS(u,t,params1),[0,sim_endtime],u(:),options);
    
    % Reshape the output array into an intelligible shape
    u = zeros([size(u) length(t)]);
    for j = 1:length(t)
        u(:,:,:,:,:,j) = reshape(u_raw(j,:),[N,N,N,3,4]);
    end
    
    % Calculate the energy in the N modes of the ROM
    energy = get_3D_energy(u,N);
    
    txt = ['\tau = ',num2str(tau)];
    plot(log(t),log(energy),'DisplayName',txt,'Color',c_list(1,:))
    box on
    xlim([log(0.1),log(sim_endtime)])
    xlabel('log(t)')
    ylabel('log(E)')
    
    % Save the results into the directory
    save(['u_n_1_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'u','-v7.3');
    save(['t_n_1_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'t');

end

hold off
legend show
saveas(gcf,sprintf('Euler_energy_%i_M_%i_N_%i_c1',sim_endtime,N_full,N),'png')


% Second order model

figure()
hold on
for i = 1:length(tau_list)
    
    tau = tau_list(i);
    
    % Find coefficients corresponding to tau and N
    ind = find(tau == tau_list); 
    c2_op_tau = c2_op{1,ind};
    ind_N = find(N == N_list); 
    c2_op_N_tau = c2_op_tau(:,ind_N);

    M = 3*N;
    
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
    params.tau = c2_op_N_tau(3);
    params2 = params;
    params2.func = @(x) t2model_RHS_tau(x);
    params2.coeff = c2_op_N_tau(1:2);
    
    % Isolate modes that do not include any components larger than N
    a = [1:N,2*M-N+1:2*M];
    
    percent_increase = 10;
    u_full = u_fullify(u,M);
    init_energy = sum(sum(sum(sum(u_full(a,a,a,:).*conj(u_full(a,a,a,:))))));
    options = odeset('AbsTol',1e-14,'Stats','on','InitialStep',1e-3,'Events',@(t,u) energy_grow(t,u,init_energy,percent_increase,params2));
    [t,u_raw] = ode45(@(t,u) RHS(u,t,params2),[0,sim_endtime],u(:),options);
    
    % Reshape the output array into an intelligible shape
    u = zeros([size(u) length(t)]);
    for j = 1:length(t)
        u(:,:,:,:,:,j) = reshape(u_raw(j,:),[N,N,N,3,4]);
    end
    
    % Calculate the energy in the N modes of the ROM
    energy = get_3D_energy(u,N);
    
    txt = ['\tau = ',num2str(tau)];
    plot(log(t),log(energy),'DisplayName',txt,'Color',c_list(1,:))
    box on
    xlim([log(0.1),log(sim_endtime)])
    xlabel('log(t)')
    ylabel('log(E)')
    
    % Save the results into the directory
    save(['u_n_2_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'u','-v7.3');
    save(['t_n_2_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'t');

end

hold off
legend show
saveas(gcf,sprintf('Euler_energy_%i_M_%i_N_%i_c2',sim_endtime,N_full,N),'png')


% Third order model

figure()
hold on
for i = 1:length(tau_list)
    
    tau = tau_list(i);
    
    % Find coefficients corresponding to tau and N
    ind = find(tau == tau_list); 
    c3_op_tau = c3_op{1,ind};
    ind_N = find(N == N_list); 
    c3_op_N_tau = c3_op_tau(:,ind_N);

    M = 3*N;
    
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
    params.tau = c3_op_N_tau(4);
    params3 = params;
    params3.func = @(x) t3model_RHS_tau(x);
    params3.coeff = c3_op_N_tau(1:3);
    
    % Isolate modes that do not include any components larger than N
    a = [1:N,2*M-N+1:2*M];
    
    percent_increase = 10;
    u_full = u_fullify(u,M);
    init_energy = sum(sum(sum(sum(u_full(a,a,a,:).*conj(u_full(a,a,a,:))))));
    options = odeset('AbsTol',1e-14,'Stats','on','InitialStep',1e-3,'Events',@(t,u) energy_grow(t,u,init_energy,percent_increase,params3));
    [t,u_raw] = ode45(@(t,u) RHS(u,t,params3),[0,sim_endtime],u(:),options);
    
    % Reshape the output array into an intelligible shape
    u = zeros([size(u) length(t)]);
    for j = 1:length(t)
        u(:,:,:,:,:,j) = reshape(u_raw(j,:),[N,N,N,3,4]);
    end
    
    % Calculate the energy in the N modes of the ROM
    energy = get_3D_energy(u,N);
    
    txt = ['\tau = ',num2str(tau)];
    plot(log(t),log(energy),'DisplayName',txt,'Color',c_list(1,:))
    box on
    xlim([log(0.1),log(sim_endtime)])
    xlabel('log(t)')
    ylabel('log(E)')
    
    % Save the results into the directory
    save(['u_n_3_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'u','-v7.3');
    save(['t_n_3_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'t');

end

hold off
legend show
saveas(gcf,sprintf('Euler_energy_%i_M_%i_N_%i_c3',sim_endtime,N_full,N),'png')


% Fourth order model

figure()
hold on
for i = 1:(length(tau_list))
    
    tau = tau_list(i);
    
    % Find coefficients corresponding to tau and N
    ind = find(tau == tau_list); 
    c4_op_tau = c4_op{1,ind};
    ind_N = find(N == N_list); 
    c4_op_N_tau = c4_op_tau(:,ind_N);

    M = 3*N;
    
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
    params.tau = c4_op_N_tau(5);
    params4 = params;
    params4.func = @(x) t4model_RHS_tau(x);
    params4.coeff = c4_op_N_tau(1:4);
    
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
    
    txt = ['\tau = ',num2str(tau)];
    plot(log(t),log(energy),'DisplayName',txt,'Color',c_list(i,:))
    box on
    xlim([log(0.1),log(sim_endtime)])
    xlabel('log(t)')
    ylabel('log(E)')
    
    % Save the results into the directory
    save(['u_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'u','-v7.3');
    save(['t_' num2str(N) '_endtime_x10_' num2str(sim_endtime*10) '_tau_x100_' num2str(tau*100) '.mat'],'t');

end

hold off
legend show
saveas(gcf,sprintf('Euler_energy_%i_M_%i_N_%i_c4',sim_endtime,N_full,N),'png')
