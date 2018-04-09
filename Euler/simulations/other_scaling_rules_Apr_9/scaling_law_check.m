function scaling_law_check(N_list,end_time,filetype)

format long
close all

load('c48_4')

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis


for i = 1:length(N_list)
  
    N = N_list(i);
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
    params.no_time = 1;
    
    
    params4 = params;
    params4.func = @(x) t4model_RHS(x);
    params4.coeff = scaling_law(N,4);
    
    if exist(sprintf('../multiple_renormalized_models_Mar_19/u_array4_%i.mat',N,end_time),'file') == 2
        
        load(sprintf('../multiple_renormalized_models_Mar_19/u_array4_%i.mat',N,end_time))
        load(sprintf('../multiple_renormalized_models_Mar_19/t4_%i',N,end_time))
        
    else
        
        % run the simulation
        options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
        [t4,u_raw4] = ode45(@(t,u) RHS(u,t,params4),[0,end_time],u(:),options);
        
        
        % reshape the output array into an intelligible shape (should make this a
        % separate function later)
        u_array4 = zeros([size(u) length(t4)]);
        for l = 1:length(t4)
            u_array4(:,:,:,:,:,l) = reshape(u_raw4(l,:),[N,N,N,3,4]);
        end
        
        save(sprintf('../multiple_renormalized_models_Mar_19/t4_%i',N,end_time),'t4');
        save(sprintf('../multiple_renormalized_models_Mar_19/u_array4_%i.mat',N,end_time),'u_array4');
        
    end
    
    
    
    params4.coeff = c4(:,i);
    
    if exist(sprintf('u_array4_%i_nofit.mat',N,end_time),'file') == 2
        
        load(sprintf('u_array4_%i_nofit.mat',N,end_time))
        load(sprintf('t4_%i_nofit',N,end_time))
        
    else
        
        % run the simulation
        options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
        [t4_nofit,u_raw4_nofit] = ode45(@(t,u) RHS(u,t,params4),[0,end_time],u(:),options);
        
        
        % reshape the output array into an intelligible shape (should make this a
        % separate function later)
        u_array4_nofit = zeros([size(u) length(t4_nofit)]);
        for l = 1:length(t4_nofit)
            u_array4_nofit(:,:,:,:,:,l) = reshape(u_raw4_nofit(l,:),[N,N,N,3,4]);
        end
        
        save(sprintf('t4_%i_nofit',N,end_time),'t4_nofit');
        save(sprintf('u_array4_%i_nofit.mat',N,end_time),'u_array4_nofit');
        
    end
    
    % plot the energy in some modes
    energy = get_3D_energy(u_array4,N);
    energy_nofit = get_3D_energy(u_array4_nofit,N);
    figure(1)
    hold off
    plot(log(t4),log(energy),'b','linewidth',2)
    hold on
    plot(log(t4_nofit),log(energy_nofit),'r','linewidth',2)
    legend(sprintf('4th order ROM for N = %i (scaling law approx)',N),sprintf('4th order ROM for N = %i (optimal coeffs)',N),'location','southwest')
    title('Energy in resolved modes','fontsize',16)
    xlabel('log(time)','fontsize',16)
    ylabel('log(energy)','fontsize',16)
    saveas(gcf,sprintf('energy_compare_%i',N),filetype)
    
    
    
    
    
    % plot the helicity
    
    w = helicity(u_array4);
    w_nofit = helicity(u_array4_nofit);
    figure(2)
    hold off
    plot(t4,w,'b','linewidth',2)
    hold on
    plot(t4_nofit,w_nofit,'r','linewidth',2)
    legend(sprintf('4th order ROM for N = %i (scaling law approx)',N),sprintf('4th order ROM for N = %i (optimal coeffs)',N),'location','southwest')
    title('Helicity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('w','fontsize',16)
    saveas(gcf,sprintf('helicity_compare_%i',N),filetype)
    
    
    
    ens = enstrophy(u_array4);
    ens_nofit = enstrophy(u_array4_nofit);
    
    figure(3)
    hold off
    plot(t4,ens,'b','linewidth',2)
    hold on
    plot(t4_nofit,ens_nofit,'r','linewidth',2)
    legend(sprintf('4th order ROM for N = %i (scaling law approx)',N),sprintf('4th order ROM for N = %i (optimal coeffs)',N),'location','northwest')
    title('Enstrophy','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('enstrophy','fontsize',16)
    saveas(gcf,sprintf('enstrophy_compare_%i',N),filetype)
    
end