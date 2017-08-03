%finds and plots scaling law behavior for t^2-model

clear all;close all;

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N_list = 32:2:60;
epsilon = fliplr(0.065:0.005:0.1);

save N_list N_list
save epsilon epsilon

t2 = zeros(length(N_list),length(epsilon));

for j = 1:length(epsilon)
    epsilon(j)
    simulation_params.epsilon = epsilon(j);  %coefficient on linear term in KdV
    simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
    simulation_params.dt = 1e-3;      %timestep
    simulation_params.endtime = 10;   %end of simulation
    simulation_params.howoften = 1;   %how often to save state vector
    simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
    simulation_params.tol = inf;    %tolerance for identifying instabilities
    simulation_params.N = 256;          %number of positive modes to simulate
    simulation_params.initial_condition = @(x) sin(x);
    
    %full model with no approximations
    simulation_params.name = 'full';
    
    [t_list,u_list] = KdV_solve(simulation_params);
    
    simulation_params = full_init(simulation_params);
    
    [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow] = generate_deriv_data_4func2(t_list,u_list,simulation_params,N_list);
    coeffs_list = no_k_dependence_coeffs2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,0);
    
    t2(:,j) = coeffs_list(:,2);
    
end

N_list_nondim = N_list*2*pi;
Re_nondim = 2*pi*(1/sqrt(2))^(1/2)./epsilon;
t2_nondim = t2./(2*sqrt(2)*pi)^2;

t2_form = find_scaling_law(t2_nondim.',Re_nondim,N_list_nondim)

save t2 t2

save t2_form t2_form



for i = 1:length(N_list)
   
    t2_eps = figure(1);
    
    hold on
    plot(log(epsilon),log(-t2(i,:)),'.','markersize',20)
    plot(log(epsilon),log(t2_form(1))+log(N_list(i))*t2_form(2)-log(epsilon)*t2_form(3),'r')
    legend('log(-t^2-coeff) from data',sprintf('log(%.3f*N^%.3f*epsilon^%.3f)',t2_form(1),t2_form(2),t2_form(3)))
    xlabel('log(epsilon)','fontsize',16)
    ylabel('log(-coeff)','fontsize',16)
    title('t^2-model epsilon scaling law','fontsize',16)
    
end


for i = 1:length(epsilon)
   
    t2_N = figure(3);
    hold on
    plot(log(N_list),log(-t2(:,i)),'.','markersize',20)
    plot(log(N_list),log(t2_form(1))+log(N_list)*t2_form(2)-log(epsilon(i))*t2_form(3),'r')
    legend('log(-t^2-coeff) from data',sprintf('log(%.3f*N^%.3f*epsilon^%.3f)',t2_form(1),t2_form(2),t2_form(3)))
    xlabel('log(N)','fontsize',16)
    ylabel('log(-coeff)','fontsize',16)
    title('t^2-model N scaling law','fontsize',16)
   
end


saveas(t2_eps,'t2_eps','png')
saveas(t2_N,'t2_N','png')
