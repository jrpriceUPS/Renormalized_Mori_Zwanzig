function [t2,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde] = t2model_term_Burgers(u_full,alpha,t0hat,t0tilde,F_modes,G_modes)
%
%Computes the complete t^2-model term of Burgers for a given state vector
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_full   =  a full state vector (positive and negative modes)
%
%  nonlin0  =  the result of a convolution of u_full with itself (markov
%              term)
%
%  uu_star  =  the unresolved part of the markov convolution
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%  G_modes  =  vector of which modes in u_full correspond to unresolved
%              modes
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t2       =  t^2-term
%
%  Ahat     =  Chat(uhat,t0hat)
%
%  Atilde   =  Ctilde(uhat,t0hat)
%
%  Bhat     =  Chat(uhat,t0tilde)
%
%  Btilde   =  Ctilde(uhat,t0tilde)
%
%  Dhat     =  Chat(t0tilde,t0tilde)
%
%  Dtilde   =  Ctilde(t0tilde,t0tilde)


A = convolution_sum_Burgers(u_full,t0hat,alpha);
Ahat = A;
Atilde = A;
Ahat(G_modes) = 0;
Atilde(F_modes) = 0;

B = convolution_sum_Burgers(u_full,t0tilde,alpha);
Bhat = B;
Btilde = B;
Bhat(G_modes) = 0;
Btilde(F_modes) = 0;

D = convolution_sum_Burgers(t0tilde,t0tilde,alpha);
Dhat = D;
Dtilde = D;
Dhat(G_modes) = 0;
Dtilde(F_modes) = 0;

term1 = convolution_sum_Burgers(u_full,Atilde-Btilde,alpha);
term1(G_modes) = 0;

t2 = 4*term1 - 2*Dhat;