function coeffs_list = k_dependent_coefficients3(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,window_width,spacing,disp_plots)
%
%Calculates optimal k-dependent 2nd order ROM coefficients for fitting both 
%the the energy in individual modes and the net energy in different sets of 
%resolved modes
%
%%%%%%%%%
%Inputs:%
%%%%%%%%%
%
%t_list  =  the list of times in the simulation we are fitting to
%
%energy_flow_list  =  an array of the exact derivatives of energy in each
%                     individual mode at each timestep
%
%nonlin0_energy_flow  =  an array of the Markov term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep
%
%nonlin1_energy_flow  =  an array of the t-model term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep
%
%nonlin2_energy_flow  =  an array of the t^2-model term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep
%
%N_list  =  the list of n resolutions for which we want to find ROM
%           coefficients
%
%window_width  =  the width of the sliding window over which we compute
%                 coefficients (set to zero if every window should begin at t=0)
%
%spacing  =  the step size in between coefficient calculations (set to
%            t_max if only want the best coefficients for a given
%            simulation)
%  
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%coeffs_list  =  an n x 2 x window array of optimal coefficients where the ith
%                column corresponds to the t^i-model

%turn off warnings for nearly singular linear solves
warning('off', 'MATLAB:nearlySingularMatrix')

%find the endtime
endtime = t_list(end);

if window_width == 0
    %if no windows are desired and only end result is wanted, set the
    %windows variable to just be [0 endtime]
    if spacing == endtime
        window_edges = endtime;
        windows = [0 endtime];
    else
        %otherwise beginning at 2, create windows of the form [0,x] for x
        %counting by spacing up to the end
        window_edges = 2:spacing:endtime;
        windows = [zeros(size(window_edges)).' window_edges.'];
    end
else
    %if sliding window is desired, construct it
    window_edges = window_width:spacing:endtime;
    windows = [window_edges.'-window_width window_edges.'];
end

%initialize output
coeffs_list = zeros(length(N_list),2*(max(N_list)-1),length(window_edges));

dt = t_list(2)-t_list(1);

for j = 1:length(window_edges)
   
    %isolate the range associated with current window
    range = round((windows(j,1)/dt:windows(j,2)/dt)+1);
    
    for i = 1:length(N_list)
        
        %solve least squares problem for fit of k-dependent coefficients
        %(minimizing error in net flow in / out of resolved modes and each 
        %individual mode)
        N = N_list(i);
        
        x = squeeze(nonlin1_energy_flow(i,2:N,range));
        y = squeeze(nonlin2_energy_flow(i,2:N,range));
        
        r = energy_flow_list(2:N,range) - squeeze(nonlin0_energy_flow(i,2:N,range));
        
        %individual modes matrix and RHS
        A1 = zeros(2*(N-1),2*(N-1));
        b1 = zeros(2*(N-1),1);
        
        for k = 1:N-1
            A1(k,1:N-1) = sum(x.*repmat(x(k,:),N-1,1),2).';
            A1(k,(N-1)+1:2*(N-1)) = sum(y.*repmat(x(k,:),N-1,1),2).';
            
            A1(k+(N-1),1:N-1) = sum(x.*repmat(y(k,:),N-1,1),2).';
            A1(k+(N-1),(N-1)+1:2*(N-1)) = sum(y.*repmat(y(k,:),N-1,1),2).';
            
            
            b1(k) = sum(sum(r,1).*x(k,:));
            b1(k+N-1) = sum(sum(r,1).*y(k,:));
        end
        
        %net energy flow matrix and RHS
        xx = sum(x.*x,2);
        xy = sum(x.*y,2);
        xr = sum(x.*r,2);
        
        yy = sum(y.*y,2);
        yr = sum(y.*r,2);
        
        
        A2 = zeros(2*(N-1),2*(N-1));
        
        A2(1:N-1,1:N-1) = diag(xx);
        A2(1:N-1,(N-1)+1:2*(N-1)) = diag(xy);
        
        A2((N-1)+1:2*(N-1),1:(N-1)) = diag(xy);
        A2((N-1)+1:2*(N-1),(N-1)+1:2*(N-1)) = diag(yy);

        
        b2 = [xr;yr];
        
        %sum to minimize both simultaneously
        A = A1 + A2;
        b = b1 + b2;
        
        coeffs = A\b;
        coeffs_list(i,1:2*(N-1),j) = coeffs;
        
        
        %produce plots if only fitting over full range
        if window_width == 0 && spacing == endtime && disp_plots == 1;
            
            %isolate coefficients
            alpha = coeffs(1:(N-1)); %t-model
            beta = coeffs((N-1)+1:2*(N-1)); %t2-model
            
            %plot coefficients as a function of wavenumber k
            figure(1)
            subplot(2,1,1)
            plot(2:N,alpha)
            ax = axis;
            axis([2,N,ax(3),ax(4)]);
            title(sprintf('N=%i t-model coefficient',N_list(i)))
            xlabel('wavenumber')
            ylabel('coefficient')
            
            subplot(2,1,2)
            plot(2:N,beta)
            ax = axis;
            axis([2,N,ax(3),ax(4)]);
            title(sprintf('N=%i t^2-model coefficient',N_list(i)))
            xlabel('wavenumber')
            ylabel('coefficient')
            
            
            
            %saveas(coeffs_subplot,sprintf('N%i_k_coeffs',N_list(i)),'png')
            
            
            
            %plot total energy flow compared against coefficient fit
            energy_plot = figure(2);
            hold off
            plot(t_list,sum(r))
            hold on
            plot(t_list,sum(repmat(alpha,1,10001).*x)+sum(repmat(beta,1,10001).*y),'r')
            title(sprintf('N = %i energy derivative',N_list(i)),'fontsize',16)
            xlabel('time','fontsize',16)
            ylabel('energy derivative','fontsize',16)
            saveas(energy_plot,sprintf('N%i_k_energy',N_list(i)),'png')
            
            %plot fit for each individual mode
            sub_modes = figure(3);
            for k = 1:N-1
                subplot(4,N/4,k)
                hold on
                plot(t_list,r(k,:))
                plot(t_list,alpha(k)*x(k,:)+beta(k)*y(k,:),'r')
            end
            saveas(sub_modes,sprintf('submodes%i',N_list(i)),'png')
        end
    end
end

