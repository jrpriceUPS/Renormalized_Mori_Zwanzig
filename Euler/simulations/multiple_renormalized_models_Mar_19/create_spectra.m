function create_spectra(N_list,skip,yes_log)

format long
close all

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

for i = 1:length(N_list)
    N = N_list(i);
    M = 3*N;
    
    load(sprintf('u_array4_%i.mat',N))
    load(sprintf('t4_%i.mat',N))
    count = 0;
    
    if ~(exist(sprintf('spect%i',N))==7)
        mkdir(sprintf('spect%i',N))
    end
    
    for j = 1:skip:length(t4)
        u_current = squeeze(u_array4(:,:,:,:,:,j));
        u_full = u_fullify(u_current,M);
        [spectrum,k_list] = energy_spectrum(u_full,M);
        
        if yes_log
            figure(1)
            hold off
            plot(log(k_list),log(spectrum),'.','markersize',20)
            title(sprintf('Energy Spectrum for N = %i at t = %i',N,t4(j)),'fontsize',16)
            xlabel('log(|k|)','fontsize',16)
            ylabel('log(energy)','fontsize',16)
            axis([0,log(N),-10,5])
            saveas(gcf,sprintf('spect%i/spectrum%i_%i',N,N,count),'png')
            count = count+1;
        else
            figure(1)
            hold off
            plot(k_list,spectrum,'.','markersize',20)
            title(sprintf('Energy Spectrum for N = %i at t = %i',N,t4(j)),'fontsize',16)
            xlabel('|k|','fontsize',16)
            ylabel('log(energy)','fontsize',16)
            axis([0,N,0,1/4])
            saveas(gcf,sprintf('spect%i/spectrum%i_%i',N,N,count),'png')
            count = count+1;
        end
    end
end