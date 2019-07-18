clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N = 24;
M = 3*N;

alpha_list = [1;0.1;0.01];
endtime_list = [2;20;200];
namelist = {'1','p1', 'p01'};
N_list = 4:2:12;

renormalization_range = zeros(2,length(alpha_list));

coefficients = zeros(length(alpha_list),4,length(N_list),4);
laws = zeros(length(alpha_list),4,2,4);
r = zeros(length(alpha_list),4,4);

coefficients_t = zeros(length(alpha_list),4,length(N_list),4);
laws_t = zeros(length(alpha_list),4,2,4);
r_t = zeros(length(alpha_list),4,4);

for j = 1:length(alpha_list)
    
    
    if ~(exist(['u_' namelist{j} '.mat'],'file') == 2)
        create_data_alpha(N,endtime_list(j),alpha_list(j),namelist{j});
    end
    %     if ~(exist(['u48_2_' namelist{j} '.mat'],'file')==2)
    %
    %         load(['u48_' namelist{j} '.mat'])
    %         load(['t48_' namelist{j} '.mat'])
    %         start_time = t(end);
    %         u0 = u(:,:,:,:,:,end);
    %
    %
    %         % make k array
    %         k_vec = [0:M-1,-M:1:-1];
    %         [kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
    %         k = zeros(2*M,2*M,2*M,3);
    %         k(:,:,:,1) = kx;
    %         k(:,:,:,2) = ky;
    %         k(:,:,:,3) = kz;
    %
    %         % load relevant parameters into parameter structure
    %         params.k = k;
    %         params.N = N;
    %         params.M = M;
    %         params.func = @(x) full_RHS(x);
    %         params.coeff = [];
    %         params.a = 2:M;
    %         params.b = 2*M:-1:M+2;
    %         params.a_tilde = N+1:M;
    %         params.a_tilde2 = 2*N+1:M;
    %         params.print_time = 1;
    %
    %         % run the simulation
    %         options = odeset('RelTol',1e-10,'Stats','on','InitialStep',1e-3);
    %         [t_new,u_raw] = ode45(@(t,u) RHS(u,t,params),[1,1.5],u0(:),options);
    %
    %         % reshape the output array into an intelligible shape (should make this a
    %         % separate function later)
    %         u_new = zeros([size(u0) length(t_new)]);
    %         for i = 1:length(t_new)
    %             u_new(:,:,:,:,:,i) = reshape(u_raw(i,:),[N,N,N,3,4]);
    %         end
    %
    %         t2 = t_new;
    %         u2 = u_new;
    %
    %         save(['u48_2_' namelist{j} '.mat'],'t2');
    %         save(['t48_2' namelist{j} '.mat'],'u2');
    %
    %     end
    
    
    
    
    load(['u_' namelist{j} '.mat'])
    load(['t_' namelist{j} '.mat'])
    
    %     load(['u48_2_' namelist{j} '.mat'])
    %     load(['t48_2' namelist{j} '.mat'])
    
    %     s = size(u);
    %
    %     u_both = zeros(s(1),s(2),s(3),s(4),s(5),length(t)+length(t2));
    %     u_both(:,:,:,:,:,1:length(t)) = u;
    %     u_both(:,:,:,:,:,length(t)+1:end) = u2;
    %
    %     t_both = [t;t2];
    
    if ~(exist(['tmodel_size_list_' namelist{j} '.mat'],'file') == 2)
        %[tmodel_size_list,tmodel_size_list_full] = resolve_array(u_both,t_both);
        [tmodel_size_list,tmodel_size_list_full] = resolve_array(u,t);
        save(['tmodel_size_list_' namelist{j} '.mat'],'tmodel_size_list')
        save(['tmodel_size_list_full_' namelist{j} '.mat'],'tmodel_size_list_full')
    end
    
    load(['tmodel_size_list_' namelist{j} '.mat'])
    min_tol = 1e-16;
    max_tol = 1e-10;
    
    viable_snapshots = find(tmodel_size_list > min_tol & tmodel_size_list < max_tol);
    
    % trim the arrays to those viable times
    u_array = u(:,:,:,:,:,viable_snapshots);
    t_array = t(viable_snapshots);
    
    % compute the coefficients for all even modes smaller than the full
    % simulation (except N = 2)
    N_list = 4:2:12;
    
    % compute the renormalization coefficients
    [c_find,laws_find,r_find] = renormalize_frac(u_array,N_list,t_array,1/4,1);
    coefficients(j,:,:,:) = c_find;
    laws(j,:,:,:) = laws_find;
    r(j,:,:) = r_find;
    save('diff_IC_coeffs_p75','coefficients')
    save('diff_IC_laws_p75','laws')
    save('diff_IC_r_p75','r')
    
    [c_find,laws_find,r_find] = renormalize_frac(u_array,N_list,t_array,3/4,1);
    coefficients(j,:,:,:) = c_find;
    laws(j,:,:,:) = laws_find;
    r(j,:,:) = r_find;
    save('diff_IC_coeffs_p25','coefficients')
    save('diff_IC_laws_p25','laws')
    save('diff_IC_r_p25','r')
    
end

markers = {'b.','r*','gs','kx'};
lines = {'b','r','g','k'};
load('diff_IC_coeffs_p75')
load('diff_IC_laws_p75')

for i = 1:length(alpha_list)
    
    for j = 1:4
        figure(j)
        subplot(2,2,1)
        plot(log(N_list),log(squeeze(coefficients(i,1,:,j))),markers{i},'markersize',20)
        hold on
        plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,1,:,j)),[1,log(N_list(end))+1]),lines{i})
        title('t-model coefficient','fontsize',16)
        xlabel('log(N)')
        ylabel('log(a_1)')
        
        if j > 1
            
            subplot(2,2,2)
            plot(log(N_list),log(squeeze(-coefficients(i,2,:,j))),markers{i},'markersize',20)
            hold on
            plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,2,:,j)),[1,log(N_list(end))+1]),lines{i})
            title('t^2-model coefficient','fontsize',16)
            xlabel('log(N)')
            ylabel('log(a_2)')
            
            if j > 2
                
                subplot(2,2,3)
                plot(log(N_list),log(squeeze(coefficients(i,3,:,j))),markers{i},'markersize',20)
                hold on
                plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,3,:,j)),[1,log(N_list(end))+1]),lines{i})
                title('t^3-model coefficient','fontsize',16)
                xlabel('log(N)')
                ylabel('log(a_3)')
                
                if j > 3
                    
                    subplot(2,2,4)
                    plot(log(N_list),log(squeeze(-coefficients(i,4,:,j))),markers{i},'markersize',20)
                    hold on
                    plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,4,:,j)),[1,log(N_list(end))+1]),lines{i})
                    title('t^4-model coefficient','fontsize',16)
                    xlabel('log(N)')
                    ylabel('log(a_4)')
                end
            end
        end
        saveas(gcf,sprintf('coeff_plot_smaller_IC_ROM_p75_%i',j),'png')
    end
end

close all




load('diff_IC_coeffs_p25')
load('diff_IC_laws_p25')

for i = 1:length(alpha_list)
    
    for j = 1:4
        figure(j)
        subplot(2,2,1)
        plot(log(N_list),log(squeeze(coefficients(i,1,:,j))),markers{i},'markersize',20)
        hold on
        plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,1,:,j)),[1,log(N_list(end))+1]),lines{i})
        title('t-model coefficient','fontsize',16)
        xlabel('log(N)')
        ylabel('log(a_1)')
        
        if j > 1
            
            subplot(2,2,2)
            plot(log(N_list),log(squeeze(-coefficients(i,2,:,j))),markers{i},'markersize',20)
            hold on
            plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,2,:,j)),[1,log(N_list(end))+1]),lines{i})
            title('t^2-model coefficient','fontsize',16)
            xlabel('log(N)')
            ylabel('log(a_2)')
            
            if j > 2
                
                subplot(2,2,3)
                plot(log(N_list),log(squeeze(coefficients(i,3,:,j))),markers{i},'markersize',20)
                hold on
                plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,3,:,j)),[1,log(N_list(end))+1]),lines{i})
                title('t^3-model coefficient','fontsize',16)
                xlabel('log(N)')
                ylabel('log(a_3)')
                
                if j > 3
                    
                    subplot(2,2,4)
                    plot(log(N_list),log(squeeze(-coefficients(i,4,:,j))),markers{i},'markersize',20)
                    hold on
                    plot([1,log(N_list(end))+1],polyval(squeeze(laws(i,4,:,j)),[1,log(N_list(end))+1]),lines{i})
                    title('t^4-model coefficient','fontsize',16)
                    xlabel('log(N)')
                    ylabel('log(a_4)')
                end
            end
        end
        saveas(gcf,sprintf('coeff_plot_smaller_IC_ROM_p25_%i',j),'png')
    end
end

close all