function coeffs_list = generate_coeffs_4figs(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N_list,window_width,spacing,disp_plots)
%
%Calculates optimal 4th order ROM coefficients for fitting the net energy
%in different sets of resolved modes
%
%%%%%%%%%
%Inputs:%
%%%%%%%%%
%
%t_list  =  the list of times in the simulation we are fitting to
%
%energy_flow_list     =  an array of the exact derivatives of energy in each
%                        individual mode at each timestep
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
%nonlin3_energy_flow  =  an array of the t^3-model term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep
%
%nonlin4_energy_flow  =  an array of the t^4-model term contribution to the
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
%disp_plots  =  a logical variable indicating whether plots should be
%               displayed during the fits (set to 1 if desired)
%  
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%coeffs_list  =  an n x 4 x window array of optimal coefficients where the ith
%                column corresponds to the t^i-model

%find the endtime
endtime = t_list(end);

if window_width == 0
    %if no windows are desired and only end result is wanted, set the
    %windows variable to just be [0 endtime]
    if spacing == endtime
        window_edges = endtime;
        windows = [0 endtime];
    else
        %otherwise beginning at 2, create windows of the form [2,x] for x
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
coeffs_list = zeros(length(N_list),4,length(window_edges));

dt = t_list(2)-t_list(1);

for j = 1:length(window_edges)
    
    %isolate the range associated with current window
    range = round((windows(j,1)/dt:windows(j,2)/dt)+1);
    
    for i = 1:length(N_list)
        
        %solve least squares problem for fit
        N = N_list(i);
        
        x = squeeze(sum(nonlin1_energy_flow(i,1:N,range)));
        y = squeeze(sum(nonlin2_energy_flow(i,1:N,range)));
        z = squeeze(sum(nonlin3_energy_flow(i,1:N,range)));
        w = squeeze(sum(nonlin4_energy_flow(i,1:N,range)));
        
        r = sum(energy_flow_list(1:N,range)).'-squeeze(sum(nonlin0_energy_flow(i,1:N,range)));
        
        A = [sum(x.*x) sum(x.*y) sum(x.*z) sum(x.*w)
            sum(y.*x) sum(y.*y) sum(y.*z) sum(y.*w)
            sum(z.*x) sum(z.*y) sum(z.*z) sum(z.*w)
            sum(w.*x) sum(w.*y) sum(w.*z) sum(w.*w)];
        
        b = [sum(x.*r);sum(y.*r);sum(z.*r);sum(w.*r)];
        
        coeffs = A\b;
        coeffs_list(i,:,j) = coeffs;
        
        %if plotting, display the current fit results
        if disp_plots == 1
            figure(1)
            subplot(3,ceil(length(N_list)/3),i)
            hold off
            plot(t_list(range),r)
            hold on
            plot(t_list(range),coeffs(1)*x + coeffs(2)*y + coeffs(3)*z + coeffs(4)*w,'r')
            current_ax = axis;
            axis([min(t_list(range)),max(t_list(range)),current_ax(3),current_ax(4)]);
            title(sprintf('N = %i energy derivative',N_list(i)))
        end
        
        
    end
end


if window_width == 0 && spacing ~= endtime
    %coeffs plot when done over the full range
    for i = 1:length(N_list)
        
        figure
        plot(window_edges,squeeze(coeffs_list(i,:,:)).')
        title(sprintf('N = %i optimal coefficients',N_list(i)),'fontsize',16)
        xlabel('Fit over range [0,x]','fontsize',16)
        ylabel('Coefficient','fontsize',16)
        saveas(gcf,sprintf('N%i_coeffs',N_list(i)),'png')
        close
        
    end
    
else
    
    %coeffs plot when done over the window
    for i = 1:length(N_list)
        
        figure
        plot(window_edges,squeeze(coeffs_list(i,:,:)).')
        title(sprintf('N = %i optimal coefficients',N_list(i)),'fontsize',16)
        xlabel(sprintf('Fit over range [x-%i,x]',window_width),'fontsize',16)
        ylabel('Coefficient','fontsize',16)
        pause
        %saveas(gcf,sprintf('N%i_coeffs',N_list(i)),'png')
        %close
        
    end
    
end