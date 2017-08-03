%a script to try out finding renormalization constants with the phase data
%included

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

clear all;close all;
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
    
    [u_deriv_list,nonlin0_deriv,nonlin1_deriv,nonlin2_deriv,nonlin3_deriv,nonlin4_deriv] = phase_deriv_data(t_list,u_list,simulation_params,N_list);
    
    %coeffs_list = phase_coeffs(t_list,u_deriv_list,nonlin0_deriv,nonlin1_deriv,nonlin2_deriv,nonlin3_deriv,nonlin4_deriv,N_list,0,10,0);
    coeffs_list = phase_coeffs2(t_list,u_deriv_list,nonlin0_deriv,nonlin2_deriv,nonlin4_deriv,N_list,0,10,0);
    
    t2(:,j) = coeffs_list(:,1);
    t4(:,j) = coeffs_list(:,2);
    
end

%compute the scaling law using log-log least squares method
t2_form = find_scaling_law(t2.',epsilon,N_list)
t4_form = find_scaling_law(t4.',epsilon,N_list)

%save t2 t2
%save t4 t4

%save t2_form t2_form
%save t4_form t4_form


%plot results
for i = 1:length(N_list)
    
    t2_eps = figure(1);
    
    hold on
    plot(log(epsilon*2^(1/4)/(2*pi)),log(-t2(i,:)/(2*sqrt(2)*pi)^2),'.','markersize',20)
    plot(log(epsilon*2^(1/4)/(2*pi)),log(t2_form(1))+log(2*pi*N_list(i))*t2_form(2)+log(epsilon*2^(1/4)/(2*pi))*t2_form(3),'r')
    legend('log(-t^2-coeff) from data',sprintf('log(%.3f*N^%.3f*epsilon^%.3f)',t2_form(1),t2_form(2),t2_form(3)))
    xlabel('log(epsilon)','fontsize',16)
    ylabel('log(-coeff)','fontsize',16)
    title('t^2-model epsilon scaling law','fontsize',16)
    
    
    t4_eps = figure(2);
    hold on
    plot(log(epsilon*2^(1/4)/(2*pi)),log(-t4(i,:)/(2*sqrt(2)*pi)^4),'.','markersize',20)
    plot(log(epsilon*2^(1/4)/(2*pi)),log(t4_form(1))+log(2*pi*N_list(i))*t4_form(2)+log(epsilon*2^(1/4)/(2*pi))*t4_form(3),'r')
    legend('log(t^4-coeff) from data',sprintf('log(%.3f*N^%.3f*epsilon^%.3f)',t4_form(1),t4_form(2),t4_form(3)))
    xlabel('log(epsilon)','fontsize',16)
    ylabel('log(-coeff)','fontsize',16)
    title('t^4-model epsilon scaling law','fontsize',16)
    
end


for i = 1:length(epsilon)
    
    t2_N = figure(3);
    hold on
    plot(log(2*pi*N_list),log(-t2(:,i)/(2*sqrt(2)*pi)^2),'.','markersize',20)
    plot(log(2*pi*N_list),log(t2_form(1))+log(2*pi*N_list)*t2_form(2)+log(epsilon(i)*2^(1/4)/(2*pi))*t2_form(3),'r')
    legend('log(-t^2-coeff) from data',sprintf('log(%.3f*N^%.3f*epsilon^%.3f)',t2_form(1),t2_form(2),t2_form(3)))
    xlabel('log(N)','fontsize',16)
    ylabel('log(-coeff)','fontsize',16)
    title('t^2-model N scaling law','fontsize',16)
    
    t4_N = figure(4);
    hold on
    plot(log(2*pi*N_list),log(-t4(:,i)/(2*sqrt(2)*pi)^4),'.','markersize',20)
    plot(log(2*pi*N_list),log(t4_form(1))+log(2*pi*N_list)*t4_form(2)+log(epsilon(i)*2^(1/4)/(2*pi))*t4_form(3),'r')
    legend('log(t^4-coeff) from data',sprintf('log(%.3f*N^%.3f*epsilon^%.3f)',t4_form(1),t4_form(2),t4_form(3)))
    xlabel('log(N)','fontsize',16)
    ylabel('log(-coeff)','fontsize',16)
    title('t^4-model N scaling law','fontsize',16)
    
end


saveas(t2_eps,'t2_eps','png')
saveas(t4_eps,'t4_eps','png')
saveas(t2_N,'t2_N','png')
saveas(t4_N,'t4_N','png')
close all



