function [nonlin0,u_full] = markov_term(u,M,N,alpha)
%
%Computes the Markov term for a given state vector and size of full model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u        =  a state vector (positive modes only)
%
%  M        =  size of the "full" model upon which we will be basing 
%              calculations
%
%  N        =  size of reduced model (I think this should be the same as
%              length(u), but I'm not 100% certain - should check out)
%
%  alpha    =  coefficient on the nonlinear term
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin0  =  the Markov term
%
%  u_full   =  the full state vector (positive and negative modes)

%fill positive and negative modes
u_full = zeros(2*M,1);
u_full(1:N) = u;
u_full(2*M-N+2:2*M) = conj(flipud(u(2:N)));

%compute first convolution
nonlin0 = convolution_sum(u_full,u_full,alpha);