function nonlin=renormalized_complete_Burgers(u,t,simulation_params)
%
%
% This function has been updated to include the renormalization exponent
% coefficient tau and clean up some of the additional parameters that are
% not needed
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
%      tau      = non-linearly of the time dependence
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin  =  nonlinear part of RHS according to this model

if simulation_params.print_time == 1
    disp(sprintf('Current time is t = %i',t))
end

% Gather parameters needed for simulation
alpha = simulation_params.alpha;
F_modes = simulation_params.F_modes;
G_modes = simulation_params.G_modes;
coeffs = simulation_params.coeffs;
N = simulation_params.N;
M = simulation_params.M;
degree = simulation_params.degree;

% Added for the ability to test different renormalization time exponents
tau = simulation_params.tau;

%compute Markov term
[t0,t0hat,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
t0 = t0(1:N);

%compute t-model term
[t1,~,~] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes);
t1 = coeffs(1)*t^(1-1*tau)*t1(1:N);

if degree > 1
    
    %compute t^2-model term
    [t2,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde] = t2model_term_Burgers(u_full,alpha,t0hat,t0tilde,F_modes,G_modes);
    t2 = coeffs(2)*t^(2-2*tau)*t2(1:N);
end


if degree > 2
    
    %compute t^3-model term
    [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde);
    t3 = coeffs(3)*t^(3-3*tau)*t3(1:N);
end


if (degree > 3) && (degree ~= 13)
    
    %compute t^4-model term
    t4 = t4model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde,Ehat,Etilde,Fhat,Ftilde);
    t4 = coeffs(4)*t^(4-4*tau)*t4(1:N);
end


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
elseif degree == 13
    nonlin = t0 + t1 + t3;
end