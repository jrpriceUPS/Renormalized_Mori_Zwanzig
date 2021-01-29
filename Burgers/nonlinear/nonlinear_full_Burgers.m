function nonlin = nonlinear_full_Burgers(u,t,simulation_params)
%
%Computes the nonlinear part of the right hand side of the full Burgers
%equation. This is a simple convolution sum.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u        =  a state vector (positive modes only)
%
%  M        =  size of system
%
%  alpha    =  coefficient on the nonlinear term
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin  =  nonlinear part of RHS of ODE

if simulation_params.print_time == 1
    disp(sprintf('Current time is t = %i',t))
end

% Gather parameters needed for simulation
alpha = simulation_params.alpha;
N = simulation_params.N;
M = simulation_params.M;

% Compute the Markov term with full model
[nonlin0,~,~,~] = markov_term_Burgers(u,M,N,alpha,[],[]);

% Retain only the positive modes we are simulating
nonlin = nonlin0(1:N);