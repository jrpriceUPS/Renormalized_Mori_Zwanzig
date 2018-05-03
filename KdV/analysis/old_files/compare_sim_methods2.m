function [timing_exact,timing_markov,timing_full,timing_final,rel_err_markov,rel_err_full,rel_err_final] = compare_sim_methods2(epsilon,N_list,fully_resolved_N,endtime,howoften)
%
%A function to compare the t^4 model against a Markov model of the same
%resolution and a Markov model that uses the same maximal FFT as the t^4
%model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  epsilon  =  degree of dispersion
%
%  N_list  =  list of ROM resolutions to check
%
%  fully_resolved_N  =  resolution needed to fully resolve system
%
%  endtime  =  end time of simulation
%
%  howoften  =  how often to save results (set to zero to save all)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  timing_exact    =  time to run exact full simulation
%
%  timing_markov   =  length(N_list) x 1 array of time to run markov model
%
%  timing_full     =  length(N_list) x 1 array of time to run markov model
%                     with same maximal FFT as ROM
%
%  timing_final    =  length(N_list) x 1 array of time to run t^4 ROM model
%
%  rel_err_markov  =  relative error of Markov model at all times
%
%  rel_err_full    =  relative error of markov model with same maximal 
%                     FFT as ROM
%
%  rel_err_final   =  relative error of t^4 ROM model


%parameter details for simulation
simulation_params.epsilon = epsilon;         %coefficient on linear term in KdV
simulation_params.alpha = 1;             %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;             %timesteps
simulation_params.endtime = endtime;          %end of simulation
simulation_params.howoften = howoften;          %how often to save state vector
simulation_params.blowup = 1;            %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = 1e10;             %tolerance for identifying instabilities
simulation_params.N = fully_resolved_N;  %number of positive modes to simulate

%full model with no approximations
model.name = 'full';
model.renormalize = 0;     %logical indicating we don't have a renormalization step here

tic;
[t_list,u_exact] = KdV_solve(simulation_params,model);
timing_exact = toc;


%generate output
timing_markov = zeros(length(N_list),1);
timing_full = zeros(length(N_list),1);
timing_final = zeros(length(N_list),1);

rel_err_markov = zeros(length(N_list),length(t_list));
rel_err_full = zeros(length(N_list),length(t_list));
rel_err_final = zeros(length(N_list),length(t_list));

%loop through all resolutions
for i = 1:length(N_list)
    time_for_loop = cputime;
    
    N = N_list(i);
    
    %Markov model of same size as maximal FFT in ROM
    model.name = 'full';
    simulation_params.N = N*2;
    
    tic
    [t_full,u_full,~] = KdV_solve(simulation_params,model);
    timing_full(i) = toc;
    
    
    %Markov model with same resolution as ROM
    simulation_params.N = N;
    
    model.name = 'full';
    
    tic
    [t_markov,u_markov,~] = KdV_solve(simulation_params,model);
    timing_markov(i) = toc;
    
    
    %ROM simulation
    model.name = 'complete_fixed';
    
    tic
    [t_final,u_final,~] = KdV_solve(simulation_params,model);
    timing_final(i) = toc;
    
    %plot energy
    figure(1)
    
    hold off
    plot(t_list,get_energy(u_exact,N_list(i)));
    hold on
    plot(t_markov,get_energy(u_markov,N_list(i)),'r');
    plot(t_full,get_energy(u_full,N_list(i)),'g');
    plot(t_final,get_energy(u_final,N_list(i)),'k');
    title(sprintf('N = %i',N_list(i)),'fontsize',16)
    legend('exact','Markov', 'Markov with same size FFT as ROM','ROM','location','northwest')
    saveas(gcf,sprintf('energy_N%i_eps%i',N_list(i),epsilon*1000),'png')
    
    %find real space results
    [~,exact] = make_real_space(u_exact(1:N_list(i),:),N_list(i));
    [~,markov] = make_real_space(u_markov(1:N_list(i),:),N_list(i));
    [~,full] = make_real_space(u_full(1:N_list(i),:),N_list(i));
    [~,final] = make_real_space(u_final(1:N_list(i),:),N_list(i));
    
    %find real space relative error
    rel_err_markov(i,1:length(t_markov)) = real(sum((markov-exact(:,1:length(t_markov))).^2)./sum(exact(:,1:length(t_markov)).^2));
    rel_err_full(i,1:length(t_full)) = real(sum((full-exact(:,1:length(t_full))).^2)./sum(exact(:,1:length(t_full)).^2));
    rel_err_final(i,1:length(t_final)) = real(sum((final-exact(:,1:length(t_final))).^2)./sum(exact(:,1:length(t_final)).^2));
    
    %plot real space relative error
    figure(2)
    hold off
    plot(t_markov,rel_err_markov(i,1:length(t_markov)),'r')
    hold on
    plot(t_full,rel_err_full(i,1:length(t_full)),'g')
    plot(t_final,rel_err_final(i,1:length(t_final)),'k')
    
    xlabel('time','fontsize',16)
    ylabel('global relative error (real space)','fontsize',16)
    title(sprintf('N = %i',N_list(i)),'fontsize',16)
    legend('Markov', 'Markov with same size FFTs as ROM','ROM','location','northwest')
    saveas(gcf,sprintf('err_N%i_eps%i',N_list(i),epsilon*1000),'png')
    
    
    figure(3)
    hold off
    plot(t_full,rel_err_full(i,1:length(t_full)),'g')
    hold on
    plot(t_final,rel_err_final(i,1:length(t_final)),'k')
    
    xlabel('time','fontsize',16)
    ylabel('global relative error (real space)','fontsize',16)
    title(sprintf('N = %i',N_list(i)),'fontsize',16)
    legend('Markov with same size FFTs as ROM','ROM','location','northwest')
    saveas(gcf,sprintf('err2_N%i_eps%i',N_list(i),epsilon*1000),'png')
    
end