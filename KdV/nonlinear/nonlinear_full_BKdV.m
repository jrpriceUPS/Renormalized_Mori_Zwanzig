function nonlin = nonlinear_full_BKdV(u,M)
%
%Computes the nonlinear part of the right hand side of the full Burgers-KdV
%equation. This is a simple convolution sum.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u        =  a state vector (positive modes only)
%
%  M        =  size of "full" model
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
[nonlin0,~] = markov_term_BKdV(u,3/2*M,M);

%retain only the positive modes we are simulating
nonlin = nonlin0(1:M);