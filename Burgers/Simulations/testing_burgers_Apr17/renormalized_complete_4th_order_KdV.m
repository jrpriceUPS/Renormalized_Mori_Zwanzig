function nonlin=renormalized_complete_4th_order_KdV(u,t,simulation_params)
%
%Computes the nonlinear part of the right hand side of the t^4-model of the
%KdV equation based upon a "full" model with M positive modes (M>N) and no
%t dependence
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
%      coeffs   =  4x1 vector of coefficients for terms in memory expansion
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
t
%gather parameters needed for simulation
alpha = simulation_params.alpha;
F_modes = simulation_params.F_modes;
G_modes = simulation_params.G_modes;
k = simulation_params.k;
coeffs = simulation_params.coeffs;
N = simulation_params.N;
M = simulation_params.M;
epsilon =simulation_params.epsilon;

degree = length(coeffs);

%compute Markov term
[nonlin0,u_full] = markov_term_KdV(u,M,N,alpha);
t0 = nonlin0(1:N);

%compute t-model term
[nonlin1,uu_star] = tmodel_term_KdV(u_full,nonlin0,alpha,F_modes);
if simulation_params.time_dependence == 1
    t1 = t*coeffs(1)*nonlin1(1:N);
else
    t1 = coeffs(1)*nonlin1(1:N);
end

if length(coeffs) > 1
    %compute t^2-model term
    [nonlin2,uk3,uu,A,A_star,B,B_star,C,C_star,D,D_star] = t2model_term_complete_KdV(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);
    if simulation_params.time_dependence == 1
        t2 = t^2*coeffs(2)*nonlin2(1:N);
    else
        t2 = coeffs(2)*nonlin2(1:N);
    end
end

if length(coeffs)>2
    
    %compute t^3-model term
    [nonlin3,uk6,E,E_star,F,F_star] = t3model_term_complete_KdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk3,A,A_star,B,B_star,C,C_star,D_star);
    if simulation_params.time_dependence == 1
        t3 = t^3*coeffs(3)*nonlin3(1:N);
    else
        t3 = coeffs(3)*nonlin3(1:N);
    end
end

if length(coeffs)>3
    
    
    %compute t^4-model term
    nonlin4 = t4model_term_complete_KdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk3,uk6,A,A_star,B,B_star,C,C_star,D,D_star,E,E_star,F,F_star);
    if simulation_params.time_dependence == 1
        t4 = t^4*coeffs(4)*nonlin4(1:N);
    else
        t4 = coeffs(4)*nonlin4(1:N);
    end
end

%compute nonlinear part of right hand side
%compute nonlinear part of right hand side
if degree == 0
    nonlin = t0;
elseif degree == 1
    nonlin = t0 + t1;
elseif degree == 2
    nonlin = t0 + t1 + t2;
elseif degree == 3
    nonlin = t0 + t1 + t2 + t3;
elseif degree == 4
    nonlin = t0 + t1 + t2 + t3 + t4;
end