function coeff_plot(N)

N_list = 2:2:N/2;
load(sprintf('c%i.mat',N))
load(sprintf('laws%i.mat',N))

for j = 1:4
    figure(j)
    subplot(2,2,1)
    plot(log(N_list),log(squeeze(coeff_array(1,:,j))),'.','markersize',20)
    hold on
    plot([1,log(N_list(end))+1],polyval(scaling_laws(1,:,j),[1,log(N_list(end))+1]),'r')
    title('t-model coefficient','fontsize',16)
    xlabel('log(N)')
    ylabel('log(a_1)')
    
    if j > 1
        
        subplot(2,2,2)
        plot(log(N_list),log(squeeze(coeff_array(2,:,j))),'.','markersize',20)
        hold on
        plot([1,log(N_list(end))+1],polyval(scaling_laws(2,:,j),[1,log(N_list(end))+1]),'r')
        title('t^2-model coefficient','fontsize',16)
        xlabel('log(N)')
        ylabel('log(a_2)')
        
        if j > 2
            
            subplot(2,2,3)
            plot(log(N_list),log(squeeze(coeff_array(3,:,j))),'.','markersize',20)
            hold on
            plot([1,log(N_list(end))+1],polyval(scaling_laws(3,:,j),[1,log(N_list(end))+1]),'r')
            title('t^3-model coefficient','fontsize',16)
            xlabel('log(N)')
            ylabel('log(a_3)')
            
            if j > 3
                
                subplot(2,2,4)
                plot(log(N_list),log(squeeze(coeff_array(4,:,j))),'.','markersize',20)
                hold on
                plot([1,log(N_list(end))+1],polyval(scaling_laws(4,:,j),[1,log(N_list(end))+1]),'r')
                title('t^4-model coefficient','fontsize',16)
                xlabel('log(N)')
                ylabel('log(a_4)')
            end
        end
    end
    saveas(gcf,sprintf('coeff_plot%i_ROM%i',N_full,j),'png')
    close
end




load(sprintf('c%i_t.mat',N))
load(sprintf('laws%i_t.mat',N))

for j = 1:4
    figure(j)
    subplot(2,2,1)
    plot(log(N_list),log(squeeze(coeff_array(1,:,j))),'.','markersize',20)
    hold on
    plot([1,log(N_list(end))+1],polyval(scaling_laws(1,:,j),[1,log(N_list(end))+1]),'r')
    title('t-model coefficient','fontsize',16)
    xlabel('log(N)')
    ylabel('log(a_1)')
    
    if j > 1
        
        subplot(2,2,2)
        plot(log(N_list),log(squeeze(coeff_array(2,:,j))),'.','markersize',20)
        hold on
        plot([1,log(N_list(end))+1],polyval(scaling_laws(2,:,j),[1,log(N_list(end))+1]),'r')
        title('t^2-model coefficient','fontsize',16)
        xlabel('log(N)')
        ylabel('log(a_2)')
        
        if j > 2
            
            subplot(2,2,3)
            plot(log(N_list),log(squeeze(coeff_array(3,:,j))),'.','markersize',20)
            hold on
            plot([1,log(N_list(end))+1],polyval(scaling_laws(3,:,j),[1,log(N_list(end))+1]),'r')
            title('t^3-model coefficient','fontsize',16)
            xlabel('log(N)')
            ylabel('log(a_3)')
            
            if j > 3
                
                subplot(2,2,4)
                plot(log(N_list),log(squeeze(coeff_array(4,:,j))),'.','markersize',20)
                hold on
                plot([1,log(N_list(end))+1],polyval(scaling_laws(4,:,j),[1,log(N_list(end))+1]),'r')
                title('t^4-model coefficient','fontsize',16)
                xlabel('log(N)')
                ylabel('log(a_4)')
            end
        end
    end
    saveas(gcf,sprintf('coeff_plot%i_ROM%i_t',N_full,j),'png')
    close
end