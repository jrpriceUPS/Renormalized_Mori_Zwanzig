function nonlin = nonlinear_full_KdV(u,M,alpha)
%
%Computes the nonlinear part of the right hand side of the full KdV
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

%compute the Markov term with full model of size 3/2 x M to dealias results
[nonlin0,~] = markov_term_KdV(u,3/2*M,M,alpha);

%retain only the positive modes we are simulating
nonlin = nonlin0(1:M);