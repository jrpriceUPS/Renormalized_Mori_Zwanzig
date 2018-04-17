function [t1,t1hat,t1tilde] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes)
%
%Computes the t-model term of KdV for a given state vector
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_full   =  a full state vector (positive and negative modes)
%
%  t0tilde  =  the result of a convolution of u_full with itself (markov
%              term)
%
%  alpha    =  coefficient on the nonlinear term
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%  G_modes  =  vector of which modes in u_full correspond to unresolved modes
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t1  =  the t-model term

%compute t model term
t1 = 2*convolution_sum_Burgers(u_full,t0tilde,alpha);
t1hat = t1;
t1tilde = t1;
t1hat(G_modes) = 0;
t1tilde(F_modes) = 0;