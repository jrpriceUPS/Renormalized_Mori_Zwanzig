function vort_int = vort_integral(N_list,end_time,filetype)

format long
close all

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

colors = linspecer(length(N_list),'qualitative');



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
    
    
        
        load(sprintf('u_array4_%i_%i.mat',N,end_time))
        load(sprintf('t4_%i_%i',N,end_time))
        
    
    
    
    
    [~,vort2] = vorticity(u_array4);
    vort_int = zeros(length(t4)-1,1);
    for j = 2:length(t4)
        vort_int(j-1) = trapz(t4(1:j),vort2(1:j).');
    end
    
    figure(5)
    hold on
    plot(t4(t4(1:end-1)<15),vort_int(t4(1:end-1)<15),'linewidth',2,'color',colors(i,:))
    legend(leg_nw{:})
    title('Integral of maximum of vorticity','fontsize',16)
    xlabel('time','fontsize',16)
    ylabel('Integral of max vorticity','fontsize',16)
    saveas(gcf,sprintf('vorticity_int_%i',N),filetype)
    
    
    
end