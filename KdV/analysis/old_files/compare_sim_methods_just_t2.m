function [timing_exact,timing_markov,timing_4,timing_2,rel_err_markov,rel_err_4,rel_err_2] = compare_sim_methods_just_t2(epsilon,N_list,fully_resolved_N,endtime,howoften)
%
%A function to compare the t^2-model against a Markov model and the
%t^4-model
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
%  timing_4        =  length(N_list) x 1 array of time to run t^4 ROM model
%
%  timing_2        =  length(N_list) x 1 array of time to run t^2 ROM model
%
%  rel_err_markov  =  relative error of Markov model at all times
%
%  rel_err_4       =  relative error of t^4 ROM model
%
%  rel_err_2       =  relative error of t^2 ROM model


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

%initialize output
timing_markov = zeros(length(N_list),1);
timing_4 = zeros(length(N_list),1);
timing_2 = zeros(length(N_list),1);

rel_err_markov = zeros(length(N_list),length(t_list));
rel_err_4 = zeros(length(N_list),length(t_list));
rel_err_2 = zeros(length(N_list),length(t_list));

%loop through resolutions
for i = 1:length(N_list)
    time_for_loop = cputime;
    
    N = N_list(i);
    
    simulation_params.N = N;
    
    %compute t^4 ROM simulation
    model.name = 'complete_fixed';
    simulation_params.order = 4;
    
    tic
    [t_4,u_4,~] = KdV_solve(simulation_params,model);
    timing_4(i) = toc;
    
    
    %compute Markov simulation
    model.name = 'full';
    
    tic
    [t_markov,u_markov,~] = KdV_solve(simulation_params,model);
    timing_markov(i) = toc;
    
    
    %compute t^2 ROM simulation
    model.name = 'complete_fixed';
    simulation_params.order = 2;
    
    tic
    [t_2,u_2,~] = KdV_solve(simulation_params,model);
    timing_2(i) = toc;
    
    %plot energy evolution
    figure(1)
    
    hold off
    plot(t_list,get_energy(u_exact,N_list(i)));
    hold on
    plot(t_markov,get_energy(u_markov,N_list(i)),'r');
    plot(t_4,get_energy(u_4,N_list(i)),'g');
    plot(t_2,get_energy(u_2,N_list(i)),'k');
    title(sprintf('N = %i',N_list(i)),'fontsize',16)
    legend('exact','Markov', '4th order ROM','2nd order ROM','location','southwest')
    saveas(gcf,sprintf('energy_N%i_eps%i',N_list(i),epsilon*1000),'png')
    
    %make real space
    [~,exact] = make_real_space(u_exact(1:N_list(i),:),N_list(i));
    [~,markov] = make_real_space(u_markov(1:N_list(i),:),N_list(i));
    [~,real_4] = make_real_space(u_4(1:N_list(i),:),N_list(i));
    [~,real_2] = make_real_space(u_2(1:N_list(i),:),N_list(i));
    
    %find real space relative error
    rel_err_markov(i,1:length(t_markov)) = real(sum((markov-exact(:,1:length(t_markov))).^2)./sum(exact(:,1:length(t_markov)).^2));
    rel_err_4(i,1:length(t_4)) = real(sum((real_4-exact(:,1:length(t_4))).^2)./sum(exact(:,1:length(t_4)).^2));
    rel_err_2(i,1:length(t_2)) = real(sum((real_2-exact(:,1:length(t_2))).^2)./sum(exact(:,1:length(t_2)).^2));
    
    %plot real space relative error
    figure(2)
    hold off
    plot(t_markov,rel_err_markov(i,1:length(t_markov)),'r')
    hold on
    plot(t_4,rel_err_4(i,1:length(t_4)),'g')
    plot(t_2,rel_err_2(i,1:length(t_2)),'k')
    
    xlabel('time','fontsize',16)
    ylabel('global relative error (real space)','fontsize',16)
    title(sprintf('N = %i',N_list(i)),'fontsize',16)
    legend('Markov', '4th order ROM','2nd order ROM','location','northwest')
    saveas(gcf,sprintf('err_N%i_eps%i',N_list(i),epsilon*1000),'png')
    
    
    figure(3)
    hold off
    plot(t_4,rel_err_4(i,1:length(t_4)),'g')
    hold on
    plot(t_2,rel_err_2(i,1:length(t_2)),'k')
    
    xlabel('time','fontsize',16)
    ylabel('global relative error (real space)','fontsize',16)
    title(sprintf('N = %i',N_list(i)),'fontsize',16)
    legend('4th order ROM','2nd order ROM','location','northwest')
    saveas(gcf,sprintf('err2_N%i_eps%i',N_list(i),epsilon*1000),'png')
    
end