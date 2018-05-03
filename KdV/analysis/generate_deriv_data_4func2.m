function [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow] = generate_deriv_data_4func2(t_list,u_list,simulation_params,N_list)
%
%Calculates rate of change of the energy exactly and due to each term in an
%ROM (up to order 2) from exact data
%
%%%%%%%%%
%Inputs:%
%%%%%%%%%
%
%t_list  =  the list of times in the simulation
%
%u_list  =  fully resolved, exact u at all times in t_list
%
%simulation_params  =  simulation parameters used to produce u_list and
%                      t_list
%
%       epsilon  =  degree of dispersion
%
%       alpha    =  coefficient of nonlinearity
%
%       N        =  number of positive resolved modes
%
%       A        =  matrix of linear part of KdV model
%
%       B        =  function handle of nonlinear part of RHS in full model
%
%N_list  =  the list of n resolutions for which we want to find ROM
%           coefficients
%  
%
%%%%%%%%%
%Output:%
%%%%%%%%%
%
%u_deriv_list  =  the exact derivative of each mode at each timestep
%
%energy_flow_list  =  an array of the exact derivatives of energy in each
%                     individual mode at each timestep
%
%nonlin0_energy_flow  =  an array of the Markov term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep
%
%nonlin1_energy_flow  =  an array of the t-model term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep
%
%nonlin2_energy_flow  =  an array of the t^2-model term contribution to the
%                        energy derivative in each individual mode at each 
%                        timestep


%load parameters for later computations
timesteps = length(t_list);
N = simulation_params.N;
alpha = simulation_params.alpha;
epsilon = simulation_params.epsilon;

%define the linear and nonlinear portions of the right hand side
A=simulation_params.A;
B=simulation_params.B;

%initialize and calculate RHS (derivative of u) for each mode at each time
u_deriv_list = zeros(N,timesteps);
for i = 1:timesteps
    u_deriv_list(:,i) = A*u_list(:,i) + B(u_list(:,i),t_list(i));
end

%compute the exact rate of change of the energy in each mode
energy_flow_list = 2*(conj(u_list).*u_deriv_list + u_list.*conj(u_deriv_list));

%initialize remaining output variables
nonlin0_energy_flow = zeros(length(N_list),N,timesteps);
nonlin1_energy_flow = zeros(length(N_list),N,timesteps);
nonlin2_energy_flow = zeros(length(N_list),N,timesteps);

%compute derivative data for each ROM and each subset of modes
for j = 1:length(N_list)
    for i = 1:timesteps
        
        %specify subset of modes and variables related to it
        N0 = N_list(j);
        u = u_list(1:N0,i);
        F_modes = [1:N0,2*N0:4*N0+2,5*N0+2:6*N0];
        G_modes = N0+1:5*N0+1;
        k = [0:3*N0-1,-3*N0:-1].';
        M = 3*N0;
        A_new = simulation_params.A(1:N0,1:N0);
        
        %compute the linear part of the rate of change of u
        linear_part = A_new*u;
        
        
        %compute Markov term
        [nonlin0,u_full] = markov_term(u,M,N0,alpha);
        
        %compute t-model term
        [nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);
        
        %compute t^2-model term
        [nonlin2,~,~,~,~,~,~,~,~,~,~] = t2model_term_complete(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);
        
        %update the Markov term to include the linear part as well
        nonlin0 = nonlin0(1:N0) + linear_part;
        
        %compute the rate of change of the energy due to each term in the
        %ROM for this subset of modes
        nonlin0_energy_flow(j,1:N0,i) = 2*(nonlin0(1:N0).*conj(u) + conj(nonlin0(1:N0)).*u);
        nonlin1_energy_flow(j,1:N0,i) = 2*(nonlin1(1:N0).*conj(u) + conj(nonlin1(1:N0)).*u);
        nonlin2_energy_flow(j,1:N0,i) = 2*(nonlin2(1:N0).*conj(u) + conj(nonlin2(1:N0)).*u);
    end
end