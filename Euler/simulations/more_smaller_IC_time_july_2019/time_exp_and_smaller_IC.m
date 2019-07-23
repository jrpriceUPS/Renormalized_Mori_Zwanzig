function time_exp_and_smaller_IC(N,endtime_list,alpha_list,scaling,coeff_plots,sim_endtime)

close all

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

markers = {'b.','r*','gs','kx'};
lines = {'b','r','g','k'};

N_list = 4:2:N/2;

coefficients = zeros(4,length(N_list),4);
laws = zeros(4,2,4);
r = zeros(4,4);

for j = 1:length(alpha_list)
    
    
    if ~(exist(['u_' num2str(alpha_list(j)) '.mat'],'file') == 2)
        create_data(N,endtime_list(j),alpha_list(j),num2str(alpha_list(j)));
    end
    
    load(['u_' num2str(alpha_list(j)) '.mat'])
    load(['t_' num2str(alpha_list(j)) '.mat'])
    
    
    
    if ~(exist(['tmodel_size_list_' num2str(alpha_list(j)) '.mat'],'file') == 2)
        [tmodel_size_list,tmodel_size_list_full] = resolve_array(u,t);
        save(['tmodel_size_list_' num2str(alpha_list(j)) '.mat'],'tmodel_size_list')
        save(['tmodel_size_list_full_' num2str(alpha_list(j)) '.mat'],'tmodel_size_list_full')
    end
    
    load(['tmodel_size_list_' num2str(alpha_list(j)) '.mat'])
    min_tol = 1e-16;
    max_tol = 1e-10;
    
    viable_snapshots = find(tmodel_size_list > min_tol & tmodel_size_list < max_tol);
    
    % trim the arrays to those viable times
    u_array = u(:,:,:,:,:,viable_snapshots);
    t_array = t(viable_snapshots);
    
    for i = 1:length(scaling)
        
        if ~(exist(['coeffs_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'],'file') == 2)
            % compute the renormalization coefficients
            [c_find,laws_find,r_find] = renormalize_frac(u_array,N_list,t_array,scaling(i),0);
            coefficients(:,:,:) = c_find;
            laws(:,:,:) = laws_find;
            r(:,:) = r_find;
            save(['coeffs_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'],'coefficients')
            save(['laws_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'],'laws')
            save(['r_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'],'r')
        end
        
        load(['coeffs_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'])
        load(['laws_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'])
        load(['r_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(N) '.mat'])
        
        if coeff_plots
            for figwin = 1:4
                figure(1)
                subplot(2,2,1)
                plot(log(N_list),log(squeeze(coefficients(1,:,figwin))),markers{figwin},'markersize',20)
                hold on
                plot([1,log(N_list(end))+1],polyval(squeeze(laws(1,:,figwin)),[1,log(N_list(end))+1]),lines{figwin})
                title(['t, eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(i)-1) 'k'],'fontsize',16)
                xlabel('log(N)')
                ylabel('log(a_1)')
                
                
                if figwin > 1
                    
                    subplot(2,2,2)
                    plot(log(N_list),log(squeeze(-coefficients(2,:,figwin))),markers{figwin},'markersize',20)
                    hold on
                    plot([1,log(N_list(end))+1],polyval(squeeze(laws(2,:,figwin)),[1,log(N_list(end))+1]),lines{figwin})
                    title(['t^2, eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(i)-1) 'k'],'fontsize',16)
                    xlabel('log(N)')
                    ylabel('log(a_2)')
                    
                    if figwin > 2
                        
                        subplot(2,2,3)
                        plot(log(N_list),log(squeeze(coefficients(3,:,figwin))),markers{figwin},'markersize',20)
                        hold on
                        plot([1,log(N_list(end))+1],polyval(squeeze(laws(3,:,figwin)),[1,log(N_list(end))+1]),lines{figwin})
                        title(['t^3, eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(i)-1) 'k'],'fontsize',16)
                        xlabel('log(N)')
                        ylabel('log(a_3)')
                        
                        if figwin > 3
                            
                            subplot(2,2,4)
                            plot(log(N_list),log(squeeze(-coefficients(4,:,figwin))),markers{figwin},'markersize',20)
                            hold on
                            plot([1,log(N_list(end))+1],polyval(squeeze(laws(4,:,figwin)),[1,log(N_list(end))+1]),lines{figwin})
                            title(['t^4, eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(i)-1) 'k'],'fontsize',16)
                            xlabel('log(N)')
                            ylabel('log(a_4)')
                        end
                    end
                end
                saveas(gcf,['coeffs_plot_eps_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(i)-1) 'k_' num2str(figwin) '.png'],'png')
            end
        end
        
        close all
    end
    
end
















%ROM simulation phase
filetype = 'png';
colors = linspecer(length(N_list),'qualitative');

ens_max = zeros(length(alpha_list),length(N_list),1);
ens_max_time = zeros(length(alpha_list),length(N_list),1);
vort_max = zeros(length(alpha_list),length(N_list),1);
vort_max_time = zeros(length(alpha_list),length(N_list),1);

bigN = N;


for l = 1:length(scaling)
    for j = 1:length(alpha_list)
        
        alpha = alpha_list(j);
        
        %initialize everything needed for the simulation
        close all
        
        % create the legends
        for i = 1:length(N_list)
            
            N = N_list(i);
            full_legend{i} = sprintf('Fourth order N = %i ROM',N);
            
        end
        
        for i = 1:length(N_list)
            
            N = N_list(i);
            
            for leg_ind = 1:i
                leg_sw{leg_ind} = full_legend{leg_ind};
                leg_se{leg_ind} = full_legend{leg_ind};
                leg_ne{leg_ind} = full_legend{leg_ind};
                leg_nw{leg_ind} = full_legend{leg_ind};
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
            eval = taylor_green(X,Y,Z)*alpha;
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
            params.time_exp = scaling(l);
            params.no_time = NaN;
            
            
            params4 = params;
            params4.func = @(x) t4model_RHS(x);
            
            zeros(length(alpha_list),4,2,4);
            
            load(['coeffs_eps_' num2str(alpha) '_scalepow_' num2str(scaling(l)-1) 'k_' num2str(bigN) '.mat'])
            params4.coeff = squeeze(coefficients(:,i,4));
            
            if exist(['u_array4_' num2str(alpha) '_scalepow_' num2str(scaling(l)-1) 'k_' num2str(N) '.mat'],'file') == 2
                
                load(['u_array4_' num2str(alpha) '_scalepow_' num2str(scaling(l)-1) 'k_' num2str(N) '.mat'])
                load(['t4_' num2str(alpha) '_scalepow_' num2str(scaling(l)-1) 'k_' num2str(N) '.mat'])
                
            else
                
                % run the simulation
                options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
                [t4,u_raw4] = ode45(@(t,u) RHS(u,t,params4),[0,sim_endtime],u(:),options);
                
                
                % reshape the output array into an intelligible shape (should make this a
                % separate function later)
                u_array4 = zeros([size(u) length(t4)]);
                for reshape_index = 1:length(t4)
                    u_array4(:,:,:,:,:,reshape_index) = reshape(u_raw4(reshape_index,:),[N,N,N,3,4]);
                end
                
                save(['t4_' num2str(alpha) '_scalepow_' num2str(scaling(l)-1) 'k_' num2str(N) '.mat'],'t4');
                save(['u_array4_' num2str(alpha) '_scalepow_' num2str(scaling(l)-1) 'k_' num2str(N) '.mat'],'u_array4');
                
            end
            
            % plot the energy in some modes
            energy = get_3D_energy(u_array4,N);
            figure(1)
            hold on
            plot(log(t4),log(energy),'linewidth',1.5,'color',colors(i,:))
            legend(leg_sw{:})
            title(['Energy , eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(l)-1) 'k'],'fontsize',16)
            xlabel('log(time)','fontsize',16)
            ylabel('log(energy)','fontsize',16)
            if sum(isnan([log(t4(2)),max(log(t4)),log(energy(end))-0.5,log(energy(1))+1]))==0
                
                if energy(1) > energy(end)
                    axis([log(t4(2)),log(sim_endtime),log(energy(end))-0.5,log(energy(1))+1])
                else
                    axis([log(t4(2)),log(sim_endtime),log(energy(1))-5,log(energy(1))+1])
                end
                
            end
            saveas(gcf,['energy_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(l)-1) 'k.' filetype],filetype)
            
            
            
            % plot the enstrophy
            ens = enstrophy(u_array4);
            [m,peak_time] = max(ens);
            ens_max(i) = m;
            ens_max_time(i) = t4(peak_time);
            
            figure(2)
            hold on
            plot(t4,ens,'linewidth',1.5,'color',colors(i,:))
            legend(leg_ne{:})
            title(['Enstrophy , eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(l)-1) 'k'],'fontsize',16)
            xlabel('time','fontsize',16)
            ylabel('enstrophy','fontsize',16)
            saveas(gcf,['enstrophy_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(l)-1) 'k.' filetype],filetype)
            
            % plot the enstrophy on a smaller domain
            figure(3)
            hold on
            plot(t4(t4<=100),ens(t4<=100),'linewidth',1.5,'color',colors(i,:))
            legend(leg_ne{:})
            title(['Enstrophy , eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(l)-1) 'k'],'fontsize',16)
            xlabel('time','fontsize',16)
            ylabel('enstrophy','fontsize',16)
            saveas(gcf,['enstrophy_trim_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(l)-1) 'k.' filetype],filetype)
            
            
            % plot the vorticity
            [~,vort2] = vorticity(u_array4);
            [m,peak_time] = max(vort2);
            vort_max(i) = m;
            vort_max_time(i) = t4(peak_time);
            
            figure(4)
            hold on
            plot(t4,vort2,'linewidth',1.5,'color',colors(i,:))
            legend(leg_ne{:})
            title(['Max vorticity, eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(l)-1) 'k'],'fontsize',16)
            xlabel('time','fontsize',16)
            ylabel('Maximal vorticity','fontsize',16)
            saveas(gcf,['vorticity_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(l)-1) 'k.' filetype],filetype)
            
            % plot the vorticity on a smaller domain
            figure(5)
            hold on
            plot(t4(t4<=100),vort2(t4<=100),'linewidth',1.5,'color',colors(i,:))
            legend(leg_ne{:})
            title(['Max vorticity, eps = ' num2str(alpha_list(j)) ', scaling = ' num2str(scaling(l)-1) 'k'],'fontsize',16)
            xlabel('time','fontsize',16)
            ylabel('maximal vorticity','fontsize',16)
            saveas(gcf,['vorticity_trim_' num2str(alpha_list(j)) '_scalepow_' num2str(scaling(l)-1) 'k.' filetype],filetype)          
        end     
    end
end
