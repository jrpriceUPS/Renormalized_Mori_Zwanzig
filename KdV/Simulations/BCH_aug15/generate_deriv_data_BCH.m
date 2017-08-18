function [u_deriv_list,energy_flow_list,nonlin0_energy_flow,nonlin1_energy_flow,nonlin2_energy_flow] = generate_deriv_data_BCH(t_list,u_list,simulation_params,N_list)
%
%Given a simulation, computes the exact energy derivative at all times.
%Also computes the contribution of each term in an order four BCH
%memory ROM to the energy.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t_list  =  vector of the times at which data were recorded
%
%  u_list  =  an Nxlength(t_list) matrix containing the solution at each
%             time in t_list (where N is the full resolution of the
%             simulation)
%
%  simulation_params: a structure containing data about the simulation
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
%  N_list  =  a list of resolutions at which to gather data
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_deriv_list  =  max_Nxtimesteps array of the exact derivative of each 
%                   mode at each timestep
%
%  energy_flow_list  =  number of resolutions x max_N x timesteps array
%                       giving the exact energy derivative of each mode at
%                       each timestep for each ROM
%
%  nonlin0_energy_flow  =  number of resolutions x max_N x timesteps array
%                          giving the R_0 contribution (Markov) to the 
%                          energy derivative of each mode at each timestep
%                          for each ROM (given the exact solution)
%
%  nonlin1_energy_flow  =  number of resolutions x max_N x timesteps array
%                          giving the R_1 contribution to the energy 
%                          derivative of each mode at each timestep for 
%                          each ROM (given the exact solution)
%
%  nonlin2_energy_flow  =  number of resolutions x max_N x timesteps array
%                          giving the R_2 contribution to the energy 
%                          derivative of each mode at each timestep for 
%                          each ROM (given the exact solution)

%find how many timesteps to take
timesteps = length(t_list);

%save maximum resolution
N = simulation_params.N;

%define the linear and nonlinear portions of the right hand side
A=simulation_params.A;
B=simulation_params.B;

%calculate the derivative of each mode at each timestep
u_deriv_list = zeros(N,timesteps);
for i = 1:timesteps
    u_deriv_list(:,i) = A*u_list(:,i) + B(u_list(:,i),t_list(i));
end

%compute the energy derivative of each mode from the mode derivatives just
%computed
energy_flow_list = 2*(conj(u_list).*u_deriv_list + u_list.*conj(u_deriv_list));


%gather parameters needed for simulation
alpha = simulation_params.alpha;
epsilon = simulation_params.epsilon;

%initialize output
nonlin0_energy_flow = zeros(length(N_list),N,timesteps);
nonlin1_energy_flow = zeros(length(N_list),N,timesteps);
nonlin2_energy_flow = zeros(length(N_list),N,timesteps);

%loop through ROM resolutions and timesteps
for j = 1:length(N_list)
    for i = 1:timesteps
        %compute each term in the memory expansion for a given ROM at a
        %given timestep using the exact solution u at that time
        
        N0 = N_list(j);
        
        u = u_list(1:N0,i);
        
        F_modes = [1:N0,2*N0:4*N0+2,5*N0+2:6*N0];
        G_modes = N0+1:5*N0+1;
        k = [0:3*N0-1,-3*N0:-1].';
        M = 3*N0;
        
        A_new = simulation_params.A(1:N0,1:N0);
        linear_part = A_new*u;
        
        
        %compute Markov term
        [nonlin0,u_full] = markov_term(u,M,N0,alpha);
        
        %compute t-model term
        [nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);
        
        %compute t^2-model term
        nonlin2 = t2model_term_BCH(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);
        
        nonlin0 = nonlin0(1:N0) + linear_part;
        
        %save energy derivative due to each term in the memory expansion
        nonlin0_energy_flow(j,1:N0,i) = 2*(nonlin0(1:N0).*conj(u) + conj(nonlin0(1:N0)).*u);
        nonlin1_energy_flow(j,1:N0,i) = 2*(nonlin1(1:N0).*conj(u) + conj(nonlin1(1:N0)).*u);
        nonlin2_energy_flow(j,1:N0,i) = 2*(nonlin2(1:N0).*conj(u) + conj(nonlin2(1:N0)).*u);
        
    end
end