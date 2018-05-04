function [slopes,slopes2,turn_times,ens_max,ens_max_time,vort_max,vort_max_time] = renormalized_multiple_res(N_list,end_time,filetype)
%
%  [slopes,slopes2,turn_times,ens_max,ens_max_time,vort_max,vort_max_time] = renormalized_multiple_res(N_list,end_time,filetype)
%
%  A function to solve Euler's equations with 4th order complete memory
%  approximation ROMs for a variety of resolutions and plot the results
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%         N  =  resolution
%
%  end_time  =  time to run simulation to
%
%  filetype  =  file type in which to save figures
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%         slopes  =  the initial log-log energy decay slopes
% 
%        slopes2  =  the secondary log-log energy decay slopes
%
%     turn_times  =  the times at which energy begins draining from the system
%
%        ens_max  =  the maximum enstrophy
%
%   ens_max_time  =  the time at which the maximum enstrophy is achieved
%
%       vort_max  =  the maximum of the vorticity
%
%  vort_max_time  =  the time at which the maximum vorticity is achieved

%initialize everything needed for the simulation
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

% create the legends
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
    if exist(sprintf('energy_%i_%i.mat',N,end_time),'file') == 2
        load(sprintf('energy_%i_%i.mat',N,end_time))
        load(sprintf('t4_%i_%i',N,end_time))
    else
        energy = get_3D_energy(u_array4,N);
        save(sprintf('energy_%i_%i.mat',N,end_time),'energy');
    end
    figure(1)
    hold on
    plot(log(t4),log(energy),'linewidth',1.5,'color',colors(i,:))
    legend(leg_sw{:})
    title('Energy in resolved modes','fontsize',16)
    xlabel('log(time)','fontsize',16)
    ylabel('log(energy)','fontsize',16)
    axis([log(t4(2)),max(log(t4)),-10,0])
    saveas(gcf,sprintf('energy_mult_%i',N),filetype)
    
    % compute slopes and energy draining times
    turn_percent = 0.1;
    start_record_percent = 0.5;
    stop_record_percent = 0.9;
    restart_percent = 0.995;
    
    energy_change = abs(energy-energy(1))/abs(energy(1));
    turn_times(i) = t4(find(energy_change>turn_percent,1));

    start_recording = t4(find(energy_change > start_record_percent,1));
    stop_recording = t4(find(energy_change > stop_record_percent,1));
    s = polyfit(log(t4(t4>start_recording&t4<stop_recording)),log(energy(t4>start_recording&t4<stop_recording)),1);
    slopes(i) = s(1);
    
    second_wind = log(t4(find(energy_change>restart_percent,1)));
    s2 = polyfit(log(t4(log(t4)>second_wind)),log(energy(log(t4)>second_wind)),1);
    slopes2(i) = s2(1);
    
    % plot the times that were singled out to ensure visually that they
    % make sense
    figure(2)
    hold on
    plot(log(t4),log(energy),'linewidth',2,'color',colors(i,:))
    plot(log(t4(find(energy_change>turn_percent,1))),log(energy(find(energy_change>turn_percent,1))),'*','markersize',20,'color',colors(i,:))
    plot(log(t4(find(energy_change>start_record_percent,1))),log(energy(find(energy_change>start_record_percent,1))),'o','markersize',20,'color',colors(i,:))
    plot(log(t4(find(energy_change>stop_record_percent,1))),log(energy(find(energy_change>stop_record_percent,1))),'o','markersize',20,'color',colors(i,:))
    plot(log(t4(find(energy_change>restart_percent,1))),log(energy(find(energy_change>restart_percent,1))),'s','markersize',20,'color',colors(i,:))
    title('Energy in resolved modes','fontsize',16)
    xlabel('log(time)','fontsize',16)
    ylabel('log(energy)','fontsize',16)
    saveas(gcf,sprintf('energy_mult_slopes_%i',N),filetype)
    
    
    % plot the enstrophy
    if exist(sprintf('enstrophy_%i_%i.mat',N,end_time),'file') == 2
        load(sprintf('enstrophy_%i_%i.mat',N,end_time))
    else
        ens = enstrophy(u_array4);
        save(sprintf('enstrophy_%i_%i.mat',N,end_time),'ens');
    end
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
    
    % plot the enstrophy on a smaller domain
    figure(4)
    hold on
    plot(t4(t4<=100),ens(t4<=100),'linewidth',1.5,'color',colors(i,:))
    legend(leg_ne{:})
    title('Enstrophy','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('enstrophy','fontsize',16)
    saveas(gcf,sprintf('enstrophy_mult_trim%i',N),filetype)
    
    
    % plot the vorticity
    if exist(sprintf('vorticity_%i_%i.mat',N,end_time),'file') == 2
        load(sprintf('vorticity_%i_%i.mat',N,end_time))
    else
        [~,vort2] = vorticity(u_array4);
        save(sprintf('vorticity_%i_%i.mat',N,end_time),'vort2');
    end
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
    
    % plot the vorticity on a smaller domain
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
    
    % plot the integral of the vorticity
    figure(7)
    hold on
    plot(t4(1:end-1),vort_int,'linewidth',1.5,'color',colors(i,:))
    legend(leg_nw{:})
    title('Integral of maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('Integral of max vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_int_%i',N),filetype)
    
    % plot the integral of the vorticity on a smaller domain
    figure(8)
    hold on
    plot(t4(t4(1:end-1)<=100),vort_int(t4(1:end-1)<=100),'linewidth',1.5,'color',colors(i,:))
    legend(leg_nw{:})
    title('Integral of maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('Integral of max vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_int_%i',N),filetype)
    
end