function coeffs_list = phase_coeffs2_only_t1_and_t2(t_list,u_deriv_list,nonlin0_deriv,nonlin2_deriv,N_list,window_width,spacing,disp_plots)
%
%Calculates optimal 2nd and 4th order ROM coefficients for fitting the net energy
%in different sets of resolved modes
%
%%%%%%%%%
%Inputs:%
%%%%%%%%%
%
%t_list  =  the list of times in the simulation we are fitting to
%
%u_deriv_list  =  an array of the exact derivatives of energy in each
%                 individual mode at each timestep
%
%nonlin0_deriv  =  an array of the Markov term contribution to the
%                  energy derivative in each individual mode at each 
%                  timestep
%
%nonlin1_deriv  =  an array of the t^1-model term contribution to the
%                  energy derivative in each individual mode at each 
%                  timestep
%
%nonlin2_deriv  =  an array of the t^2-model term contribution to the
%                  energy derivative in each individual mode at each 
%                  timestep
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
%coeffs_list  =  an n x 2 x window array of optimal coefficients where the ith
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
coeffs_list = zeros(length(N_list),length(window_edges));

dt = t_list(2)-t_list(1);

for j = 1:length(window_edges)
    
    %isolate the range associated with current window
    range = round((windows(j,1)/dt:windows(j,2)/dt)+1);
    
    for i = 1:length(N_list)
        
        %solve least squares problem for fit
        N = N_list(i);
        
        w_complex = squeeze(nonlin2_deriv(i,1:N,range));
        
        r_complex = u_deriv_list(1:N,range)-squeeze(nonlin0_deriv(i,1:N,range));
        
        w = real(w_complex);
        r = real(r_complex);
        
        A1 = sum(sum(w.*w));
          
        b1 = sum(sum(w.*r));
        
        w = imag(w_complex);
        r = imag(r_complex);
        
        A2 = sum(sum(w.*w));
          
        b2 = sum(sum(w.*r));
        
        A = A1 + A2;
        b = b1 + b2;
        
        coeffs = A\b;
        coeffs_list(i,j) = coeffs;
        
        %if plotting, display the current fit results
        if disp_plots == 1
            figure(1)
            subplot(2,1,1)
            hold off
            plot(t_list(range),real(sum(r_complex,1)))
            hold on
            plot(t_list(range),sum(coeffs*real(w_complex),1),'r')
            current_ax = axis;
            axis([min(t_list(range)),max(t_list(range)),current_ax(3),current_ax(4)]);
            subplot(2,1,2)
            hold off
            plot(t_list(range),imag(sum(r_complex,1)))
            hold on
            plot(t_list(range),sum(coeffs*imag(w_complex),1),'r')
            current_ax = axis;
            axis([min(t_list(range)),max(t_list(range)),current_ax(3),current_ax(4)]);
            title(sprintf('N = %i net energy derivative',N_list(i)))
            legend('exact','ROM fit')
            saveas(gcf,sprintf('total_energy_%i',N),'png')
            
            figure(2)
            for k = 1:N-1
                subplot(4,ceil(N/4),k)
                hold off
                plot(t_list(range),real(r_complex(k+1,:)))
                hold on
                plot(t_list(range),real(coeffs*w_complex(k+1,:)),'r')
                title(sprintf('N = %i',k+1))
            end
            
            figure(3)
            for k = 1:N-1
                subplot(4,ceil(N/4),k)
                hold off
                plot(t_list(range),imag(r_complex(k+1,:)))
                hold on
                plot(t_list(range),imag(coeffs*w_complex(k+1,:)),'r')
                title(sprintf('N = %i',k+1))
            end
      
            saveas(gcf,sprintf('individual_modes_%i',N),'png')
        end
        
        
    end
end