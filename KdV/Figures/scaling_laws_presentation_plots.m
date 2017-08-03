%produces plots demonstrating impressive scaling law behavior

%produces coefficients, finds scaling law behavior for them, and plots the
%results to demonstrate how well they fit!

clear all;close all;

addpath ../simulation_functions
addpath ../nonlinear
addpath ../analysis

%use N values that will yield valid ROMs for all chosen epsilon values
N_list = 32:2:60;
epsilon = fliplr(0.065:0.005:0.1);

save N_list N_list
save epsilon epsilon

t2 = zeros(length(N_list),length(epsilon));
t4 = zeros(length(N_list),length(epsilon));

%run exact solution to time 10 and use to find t2 and t4 coefficients for
%each ROM
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
    
    [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow] = generate_deriv_data_4func(t_list,u_list,simulation_params,N_list);
    coeffs_list = no_k_dependence_coeffs(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,nonlin3_energy_flow,nonlin4_energy_flow,N_list,0,10,0);
    
    t2(:,j) = coeffs_list(:,2);
    t4(:,j) = coeffs_list(:,4);
    
end


t2_form = find_scaling_law(t2.',epsilon,N_list);
t4_form = find_scaling_law(t4.',epsilon,N_list);


%plot scaling laws
results_style2 = {'.','o','*','x','+'};
results_style4 = {'k.','ko','k*','kx','k+'};

%plot results
j = 0;
for i = 1:3:length(N_list)
    j = j+1;
   
    t_eps = figure(1);
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    
    semilogy(log(epsilon),-t2(i,:),results_style2{j},'markersize',10)
    hold on
    loglog(log(epsilon),t2_form(1)*N_list(i)^t2_form(2)*epsilon.^t2_form(3),'b','linewidth',2)
    semilogy(log(epsilon),-t4(i,:),results_style4{j},'markersize',10)
    loglog(log(epsilon),t4_form(1)*N_list(i)^t4_form(2)*epsilon.^t4_form(3),'k','linewidth',2)
    xlabel('log(epsilon)')
    ylabel('log(|alpha_i|)')
    
end
ax = axis;
axis([-2.8,-2.2,ax(3),10*ax(4)])

j = 0;
for i = 1:2:length(epsilon)
    j = j+1;
    
    t_N = figure(2);
    set(gca,'FontSize',16)
    set(gca,'FontSize',16)
    semilogy(log(N_list),-t2(:,i),results_style2{j},'markersize',10)
    xlabel('log(N)')
    ylabel('log(|alpha_i|)')
    
    hold on
    loglog(log(N_list),t2_form(1)*N_list.^t2_form(2)*epsilon(i)^t2_form(3),'b','linewidth',2)
    semilogy(log(N_list),-t4(:,i),results_style4{j},'markersize',10)
    loglog(log(N_list),t4_form(1)*N_list.^t4_form(2)*epsilon(i)^t4_form(3),'k','linewidth',2)
    
end
ax = axis;
axis([3.4,4.2,ax(3),10*ax(4)])

saveas(t_eps,'t_eps_pres','png')
saveas(t_N,'t_N_pres','png')
close all