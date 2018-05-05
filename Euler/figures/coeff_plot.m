function coeff_plot(N)
%
% coeff_plot(N)
%
% Computes constant and algebraically decaying renormalization coefficients
% for Euler's equations and plot them
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%  N  =  resolution of exact model to use

format long
close all

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

% load exact data and use 4:2:N/2 for the ROM resolutions to fit
N_list = 4:2:N/2;
load(sprintf('c%i.mat',N))
load(sprintf('laws%i.mat',N))

coeff_array = c48;
scaling_laws = laws48;

% calculate and plot the optimal scaling laws and plot the results
colors = linspecer(4,'qualitative');
dots = zeros(4,4);

for j = 1:4
    figure(1)
    subplot(2,2,1)
    a = plot(log(N_list),log(squeeze(coeff_array(1,:,j))),'.','markersize',20,'color',colors(j,:));
    dots(1,j) = a;
    hold on
    plot([1,log(N_list(end))*1.1],polyval(scaling_laws(1,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
    
    title('t-model coefficient','fontsize',16)
    xlabel('log(N)')
    ylabel('log(a_1)')
    
    
    if j > 1
        
        subplot(2,2,2)
        a = plot(log(N_list),log(squeeze(-coeff_array(2,:,j))),'.','markersize',20,'color',colors(j,:));
        dots(2,j) = a;
        hold on
        plot([1,log(N_list(end))*1.1],polyval(scaling_laws(2,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
        
        title('t^2-model coefficient','fontsize',16)
        xlabel('log(N)')
        ylabel('log(a_2)')
        
        if j > 2
            
            subplot(2,2,3)
            a = plot(log(N_list),log(squeeze(coeff_array(3,:,j))),'.','markersize',20,'color',colors(j,:));
            dots(3,j) = a;
            hold on
            plot([1,log(N_list(end))*1.1],polyval(scaling_laws(3,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
            
            title('t^3-model coefficient','fontsize',16)
            xlabel('log(N)')
            ylabel('log(a_3)')
            
            if j > 3
                
                subplot(2,2,4)
                a = plot(log(N_list),log(squeeze(-coeff_array(4,:,j))),'.','markersize',20,'color',colors(j,:));
                dots(4,j) = a;
                hold on
                plot([1,log(N_list(end))*1.1],polyval(scaling_laws(4,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
                
                title('t^4-model coefficient','fontsize',16)
                xlabel('log(N)')
                ylabel('log(a_4)')
            end
        end
    end
end

subplot(2,2,1)
legend([dots(1,1) dots(1,2) dots(1,3) dots(1,4)],{'n = 1','n = 2','n = 3','n = 4'},'location','southwest')

subplot(2,2,2)
legend([dots(2,2) dots(2,3) dots(2,4)],{'n = 2','n = 3','n = 4'},'location','southwest')

subplot(2,2,3)
legend([dots(3,3) dots(3,4)],{'n = 3','n = 4'},'location','southwest')

subplot(2,2,4)
legend([dots(4,4)],{'n = 4'},'location','southwest')

saveas(gcf,sprintf('coeff_plot%i_ROM',N),'png')
close



% load the exact data
load(sprintf('c%i_t.mat',N))
load(sprintf('laws%i_t.mat',N))

coeff_array = c48_t;
scaling_laws = laws48_t;

for j = 1:4
    figure(1)
    subplot(2,2,1)
    a = plot(log(N_list),log(squeeze(coeff_array(1,:,j))),'.','markersize',20,'color',colors(j,:));
    dots(1,j) = a;
    hold on
    plot([1,log(N_list(end))*1.1],polyval(scaling_laws(1,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
    
    title('t-model coefficient','fontsize',16)
    xlabel('log(N)')
    ylabel('log(a''_1)')
    
    if j > 1
        
        subplot(2,2,2)
        a = plot(log(N_list),log(squeeze(-coeff_array(2,:,j))),'.','markersize',20,'color',colors(j,:));
        dots(2,j) = a;
        hold on
        plot([1,log(N_list(end))*1.1],polyval(scaling_laws(2,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
        
        title('t^2-model coefficient','fontsize',16)
        xlabel('log(N)')
        ylabel('log(a''_2)')
        
        if j > 2
            
            subplot(2,2,3)
            a = plot(log(N_list),log(squeeze(coeff_array(3,:,j))),'.','markersize',20,'color',colors(j,:));
            dots(3,j) = a;
            hold on
            plot([1,log(N_list(end))*1.1],polyval(scaling_laws(3,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
            
            title('t^3-model coefficient','fontsize',16)
            xlabel('log(N)')
            ylabel('log(a''_3)')
            
            if j > 3
                
                subplot(2,2,4)
                a = plot(log(N_list),log(squeeze(-coeff_array(4,:,j))),'.','markersize',20,'color',colors(j,:));
                dots(4,j) = a;
                hold on
                plot([1,log(N_list(end))*1.1],polyval(scaling_laws(4,:,j),[1,log(N_list(end))*1.1]),'color',colors(j,:))
                
                title('t^4-model coefficient','fontsize',16)
                xlabel('log(N)')
                ylabel('log(a''_4)')
            end
        end
    end
end

% calculate and plot the optimal scaling laws and plot the results
subplot(2,2,1)
legend([dots(1,1) dots(1,2) dots(1,3) dots(1,4)],{'n = 1','n = 2','n = 3','n = 4'},'location','southwest')

subplot(2,2,2)
legend([dots(2,2) dots(2,3) dots(2,4)],{'n = 2','n = 3','n = 4'},'location','southwest')

subplot(2,2,3)
legend([dots(3,3) dots(3,4)],{'n = 3','n = 4'},'location','southwest')

subplot(2,2,4)
legend([dots(4,4)],{'n = 4'},'location','southwest')
saveas(gcf,sprintf('coeff_plot%i_ROM_t',N),'png')
close