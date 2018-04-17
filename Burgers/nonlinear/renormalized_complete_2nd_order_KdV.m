function nonlin=renormalized_complete_2nd_order_KdV(u,t,simulation_params)
%
%Computes the nonlinear part of the right hand side of the complete 
%t^2-model of the KdV equation based upon a "full" model with M positive 
%modes (M>N) with no time dependence and fixed coefficients
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  current state vector
%
%  t  =  current time
%
%  simulation_params: structure containing details of simulation
%
%      alpha    =  coefficient on the nonlinear term
%
%      F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%      G_modes  =  vector of which modes in u_full correspond to unresolved
%                  modes
%
%      k        =  vector of wavenumbers corresponding to entries of u_full
%
%      coeffs   =  2x1 vector of coefficients for terms in memory expansion
%
%      N        =  resolution of ROM
%
%      M        =  resolution of "full" model
%
%      epsilon  =  degree of dispersion
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin  =  nonlinear part of RHS according to this model



%gather parameters needed for simulation
alpha = simulation_params.alpha;
F_modes = simulation_params.F_modes;
G_modes = simulation_params.G_modes;
k = simulation_params.k;
coeffs = simulation_params.coeffs;
N = simulation_params.N;
M = simulation_params.M;
epsilon =simulation_params.epsilon;

%compute Markov term
    [nonlin0,u_full] = markov_term_KdV(u,M,N,alpha);

%compute t-model term
    [nonlin1,uu_star] = tmodel_term_KdV(u_full,nonlin0,alpha,F_modes);

%compute t^2-model term
    [nonlin2,~,~,~,~,~,~,~,~,~] = t2model_term_complete_KdV(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute nonlinear part of right hand side of t-model for KdV
if simulation_params.time_dependence == 1
    
    nonlin = nonlin0(1:N) + t*coeffs(1)*nonlin1(1:N) - t^2/2*coeffs(2)*nonlin2(1:N);
    
elseif simulation_params.time_dependence == 0
    
    nonlin = nonlin0(1:N) + coeffs(1)*nonlin1(1:N) + coeffs(2)*nonlin2(1:N);
    
end