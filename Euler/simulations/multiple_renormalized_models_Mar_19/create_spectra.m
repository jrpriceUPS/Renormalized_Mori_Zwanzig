function create_spectra(N_list)

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
    
    for j = 1:length(t4)
        u_current = squeeze(u_array4(:,:,:,:,:,j));
        u_full = u_fullify(u_current,M);
        [spectrum,k_list] = energy_spectrum(u_full,M);
        figure(1)
        hold off
        plot(log(k_list),log(spectrum),'.','markersize',20)
        title(sprintf('Energy Spectrum for N = %i at t = %i',N,t4(j)),'fontsize',16)
        xlabel('log(|k|^2)','fontsize',16)
        ylabel('log(energy)','fontsize',16)
        saveas(gcf,sprintf('spectrum%i_%i',N,j),'png')
    end
end