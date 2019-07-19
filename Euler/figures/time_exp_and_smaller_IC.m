function time_exp_and_smaller_IC(N,endtime_list,alpha_list,scaling,coeff_plots)

close all

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

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

