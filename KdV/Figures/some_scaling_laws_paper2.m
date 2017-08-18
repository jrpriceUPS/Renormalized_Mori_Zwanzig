%produces coefficients for only up to t2, finds scaling law behavior

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

%run exact solution to time 10 and use to find t2 coefficients for
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
    coeffs_list = no_k_dependence_coeffs2(t_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow,N_list,0,10,0);
    
    t2(:,j) = coeffs_list(:,2);
    
end

%compute the scaling law using log-log least squares method
N_list_nondim = N_list*2*pi;
Re_nondim = 2*pi*(1/sqrt(2))^(1/2)./epsilon;
t2_nondim = t2./(2*sqrt(2)*pi)^2;

t2_form = find_scaling_law(t2_nondim.',Re_nondim,N_list_nondim)


results_style2 = {'.','o','*','x','+'};

%plot results
j = 0;
for i = 1:3:length(N_list)
    j = j+1;
    
    stacked_plots = figure(1);
    subplot(1,2,1)
    set(gca,'FontSize',16)
    
    if j == 1
        semilogy(log(Re_nondim),-t2_nondim(i,:),results_style2{j},'markersize',30)
    else
        semilogy(log(Re_nondim),-t2_nondim(i,:),results_style2{j},'markersize',10)
    end
    hold on
    loglog(log(Re_nondim),t2_form(1)*N_list_nondim(i)^t2_form(2)*Re_nondim.^t2_form(3),'b')
    
    xlabel('log(Re)')
    ylabel('log(|Pi_i|)')
    
end
ax = axis;
axis([3.9,4.5,ax(3),10*ax(4)])

j = 0;
for i = 1:2:length(Re_nondim)
    j = j+1;
    
    subplot(1,2,2)
    set(gca,'FontSize',16)
    if j == 1
        semilogy(log(N_list_nondim),-t2_nondim(:,i),results_style2{j},'markersize',30)
    else
        semilogy(log(N_list_nondim),-t2_nondim(:,i),results_style2{j},'markersize',10)
    end
    xlabel('log(Lambda)')
    ylabel('log(|Pi_i|)')
    
    hold on
    loglog(log(N_list_nondim),t2_form(1)*N_list_nondim.^t2_form(2)*Re_nondim(i)^t2_form(3),'b')
    
    
end
ax = axis;
axis([5.2,6,ax(3),10*ax(4)])