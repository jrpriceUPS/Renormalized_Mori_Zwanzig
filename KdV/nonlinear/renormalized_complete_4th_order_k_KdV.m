function nonlin=renormalized_complete_4th_order_k_KdV(u,t,simulation_params)
%
%Computes the nonlinear part of the right hand side of the t^4-model of the
%KdV equation based upon a "full" model with M positive modes (M>N) with no
%time-dependence in the coefficients, and different coefficients for each
%mode (will also do t^2-model if coeff vector is different size
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
%      coeffs   =  4(N-1) x 1 vector of coefficients for terms in memory expansion
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
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);

%compute t^2-model term
[nonlin2,uk3,uu,A,A_star,B,B_star,C,C_star,D,D_star] = t2model_term_complete(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

if length(coeffs) == 4*(N-1)
    %compute t^3-model term
    [nonlin3,uk6,E,E_star,F,F_star] = t3model_term_complete(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk3,A,A_star,B,B_star,C,C_star,D_star);
end

if length(coeffs) == 4*(N-1)
    %compute t^4-model term
    nonlin4 = t4model_term_complete(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk3,uk6,A,A_star,B,B_star,C,C_star,D,D_star,E,E_star,F,F_star);
end


%compute nonlinear part of right hand side
if length(coeffs) == 4*(N-1)
    nonlin = nonlin0(1:N) + [0;coeffs(1:N-1).'].*nonlin1(1:N) + [0;coeffs((N-1)+1:2*(N-1)).'].*nonlin2(1:N) + [0;coeffs(2*(N-1)+1:3*(N-1)).'].*nonlin3(1:N) + [0;coeffs(3*(N-1)+1:4*(N-1)).'].*nonlin4(1:N);
else
    nonlin = nonlin0(1:N) + [0;coeffs(1:N-1).'].*nonlin1(1:N) + [0;coeffs((N-1)+1:2*(N-1)).'].*nonlin2(1:N);
end