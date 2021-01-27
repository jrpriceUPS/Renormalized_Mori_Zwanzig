function [t0,t0hat,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes)
%
%Computes the Markov term of Burgers for a given state vector and size of full 
%model
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
%  F_modes  =  a cell array of indices for resolved modes in the full
%              model
%
%  G_modes  =  a cell array of indices for unresolved modes in the full
%              model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t0       =  the Markov term
% 
%  t0hat    =  resolved modes of Markov term
%
%  t0tilde  =  unresolved modes of Markov term
%
%  u_full   =  the full state vector (positive and negative modes)

%fill positive and negative modes - set -N mode to zero
u_full = zeros(2*M,1);
u_full(1:N) = u;
u_full(2*M-N+2:2*M) = conj(flipud(u(2:N)));

%compute first convolution
t0 = convolution_sum_Burgers(u_full,u_full,alpha);
t0hat = t0;
t0tilde = t0;
t0hat(G_modes) = 0;
t0tilde(F_modes) = 0;