function [slopes,slopes2,turn_times,ens_max,ens_max_time,vort_max,vort_max_time] = renormalized_multiple_res(N_list,end_time,filetype)

format long
close all

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

colors = linspecer(length(N_list),'qualitative');

slopes = zeros(length(N_list),1);
slopes2 = zeros(length(N_list),1);
turn_times = zeros(length(N_list),1);
ens_max = zeros(length(N_list),1);
ens_max_time = zeros(length(N_list),1);
vort_max = zeros(length(N_list),1);
vort_max_time = zeros(length(N_list),1);

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
        leg_nw{j} = full_legend{j};
    end
    
    leg_sw{i+1} = 'location';
    leg_sw{i+2} = 'southwest';
    leg_se{i+1} = 'location';
    leg_se{i+2} = 'southeast';
    leg_ne{i+1} = 'location';
    leg_ne{i+2} = 'northeast';
    leg_nw{i+1} = 'location';
    leg_nw{i+2} = 'northwest';
    
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
    save(sprintf('energy_%i_%i',N,end_time),'energy')
    figure(1)
    hold on
    plot(log(t4),log(energy),'linewidth',1.5,'color',colors(i,:))
    legend(leg_sw{:})
    title('Energy in resolved modes','fontsize',16)
    xlabel('log(time)','fontsize',16)
    ylabel('log(energy)','fontsize',16)
    axis([min(log(t4)),max(log(t4)),energy(end),0])
    saveas(gcf,sprintf('energy_mult_%i',N),filetype)
    
    energy_change = abs(energy-energy(1))/abs(energy(1));
    turn_times(i) = t4(find(energy_change>0.10,1));

    
    stop_recording = log(t4(find(energy_change>0.90,1)));
    s = polyfit(log(t4(log(t4)>turn_times(i)&log(t4)<stop_recording)),log(energy(log(t4)>turn_times(i)&log(t4)<stop_recording)),1);
    slopes(i) = s(1);
    
    second_wind = log(t4(find(energy_change>0.99,1)));
    s2 = polyfit(log(t4(log(t4)>second_wind)),log(energy(log(t4)>second_wind)),1);
    slopes2(i) = s2(1);
    
    figure(2)
    hold on
    plot(log(t4),log(energy),'linewidth',2,'color',colors(i,:))
    plot(log(t4(find(energy_change>0.10,1))),log(energy(find(energy_change>0.10,1))),'o','markersize',20,'color',colors(i,:))
    plot(log(t4(find(energy_change>0.90,1))),log(energy(find(energy_change>0.90,1))),'o','markersize',20,'color',colors(i,:))
    plot(log(t4(find(energy_change>0.99,1))),log(energy(find(energy_change>0.99,1))),'s','markersize',20,'color',colors(i,:))
    title('Energy in resolved modes','fontsize',16)
    xlabel('log(time)','fontsize',16)
    ylabel('log(energy)','fontsize',16)
    saveas(gcf,sprintf('energy_mult_slopes_%i',N),filetype)
    
    
    
    ens = enstrophy(u_array4);
    [m,peak_time] = max(ens);
    ens_max(i) = m;
    ens_max_time(i) = t4(peak_time);
    
    figure(3)
    hold on
    plot(t4,ens,'linewidth',1.5,'color',colors(i,:))
    legend(leg_ne{:})
    title('Enstrophy','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('enstrophy','fontsize',16)
    saveas(gcf,sprintf('enstrophy_mult_%i',N),filetype)
    
    figure(4)
    hold on
    plot(t4(t4<=100),ens(t4<=100),'linewidth',1.5,'color',colors(i,:))
    legend(leg_ne{:})
    title('Enstrophy','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('enstrophy','fontsize',16)
    saveas(gcf,sprintf('enstrophy_mult_trim%i',N),filetype)
    
    [~,vort2] = vorticity(u_array4);
    [m,peak_time] = max(vort2);
    vort_max(i) = m;
    vort_max_time(i) = t4(peak_time);
    
    figure(5)
    hold on
    plot(t4,vort2,'linewidth',1.5,'color',colors(i,:))
    legend(leg_ne{:})
    title('Maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('Maximal vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_mult_%i',N),filetype)
    
    figure(6)
    hold on
    plot(t4(t4<=100),vort2(t4<=100),'linewidth',1.5,'color',colors(i,:))
    legend(leg_ne{:})
    title('Maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('maximal vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_mult_trim%i',N),filetype)
    
    vort_int = zeros(length(t4)-1,1);
    for j = 2:length(t4)
        vort_int(j-1) = trapz(t4(1:j),vort2(1:j).');
    end
    
    figure(7)
    hold on
    plot(t4(1:end-1),vort_int,'linewidth',1.5,'color',colors(i,:))
    legend(leg_nw{:})
    title('Integral of maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('Integral of max vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_int_%i',N),filetype)
    
    figure(8)
    hold on
    plot(t4(t4(1:end-1)<=100),vort_int(t4(1:end-1)<=100),'linewidth',1.5,'color',colors(i,:))
    legend(leg_nw{:})
    title('Integral of maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('Integral of max vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_int_%i',N),filetype)
    
end