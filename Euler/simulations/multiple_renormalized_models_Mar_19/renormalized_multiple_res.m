function [slopes,turn_times] = renormalized_multiple_res(N_list,end_time,filetype)

format long
close all

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

colors = linspecer(length(N_list),'qualitative');

slopes = zeros(length(N_list),1);
turn_times = zeros(length(N_list),1);

for i = 1:length(N_list)
    
    N = N_list(i);
    full_legend{i} = sprintf('Fourth order N = %i ROM',N);
    
end

for i = 1:length(N_list)
    
    N = N_list(i);
    
    for j = 1:i
        leg_sw{j} = full_legend{j};
        leg_se{j} = full_legend{j};
        leg_ne{j} = full_legend{j};
    end
    
    leg_sw{i+1} = 'location';
    leg_sw{i+2} = 'southwest';
    leg_se{i+1} = 'location';
    leg_se{i+2} = 'southeast';
    leg_ne{i+1} = 'location';
    leg_ne{i+2} = 'northeast';
    
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
    
    if exist(sprintf('u_array4_%i_%i.mat',N,end_time),'file') == 2
        
        load(sprintf('u_array4_%i_%i.mat',N,end_time))
        load(sprintf('t4_%i_%i',N,end_time))
        
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
        
        save(sprintf('t4_%i_%i',N,end_time),'t4');
        save(sprintf('u_array4_%i_%i.mat',N,end_time),'u_array4');
        
    end
    
    % plot the energy in some modes
    energy = get_3D_energy(u_array4,N);
    figure(1)
    hold on
    plot(log(t4),log(energy),'linewidth',2,'color',colors(i,:))
    legend(leg_sw{:})
    title('Energy in resolved modes','fontsize',16)
    xlabel('log(time)','fontsize',16)
    ylabel('log(energy)','fontsize',16)
    saveas(gcf,sprintf('energy_mult_%i',N),filetype)
    
    energy_change = abs(log(energy)-log(energy(1)))/abs(log(energy(1)));
    turn_times(i) = log(t4(find(energy_change>1e-1,1)));
    current_turn_times = turn_times(1:i)
    save current_turn_times current_turn_times
    
    s = polyfit(log(t4(log(t4)>1.1*turn_times(i))),log(energy(log(t4)>1.1*turn_times(i))),1);
    slopes(i) = s(1);
    
    
    
    current_slopes = slopes(1:i)
    save current_slopes current_slopes
    
    
    
    
    
    % plot the helicity
    
    w = helicity(u_array4);
    figure(2)
    hold on
    plot(t4,w,'linewidth',2,'color',colors(i,:))
    legend(leg_sw{:})
    title('Helicity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('w','fontsize',16)
    saveas(gcf,sprintf('helicity_mult_%i',N),filetype)
    
    
%     if exist(sprintf('d4_%i_%i.mat',N,end_time),'file') == 2
%         load(sprintf('d4_%i_%i',N,end_time))
%         d = d4;
%     else
%         d = energy_derivative(u_array4,t4,params4);
%         d4 = d;
%         save(sprintf('d4_%i_%i',N,end_time),'d4')
%     end
%     figure(3)
%     hold on
%     plot(t4,d,'linewidth',2,'color',colors(i,:))
%     legend(leg_se{:})
%     title('Energy Derivative','fontsize',16)
%     xlabel('time','fontsize',16)
%     ylabel('w','fontsize',16)
%     saveas(gcf,sprintf('energy_deriv_mult_%i',end_time),filetype)
    
    
    
    ens = enstrophy(u_array4);
    
    figure(4)
    hold on
    plot(t4,ens,'linewidth',2,'color',colors(i,:))
    legend(leg_ne{:})
    title('Enstrophy','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('enstrophy','fontsize',16)
    saveas(gcf,sprintf('enstrophy_mult_%i',N),filetype)
    
end