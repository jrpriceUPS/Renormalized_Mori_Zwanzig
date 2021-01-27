function [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde)
%
%Computes the complete t^3-model term of Burgers for a given state vector
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  alpha    =  coefficient on the nonlinear term
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%  G_modes  =  vector of which modes in u_full correspond to unresolved
%              modes
%
%  u_full   =  a full state vector (positive and negative modes)
%
%  t0hat    =  resolved part of Markov term
%
%  t0tilde  =  unresolved part of Markov term
%
%  Ahat     =  Chat(uhat,t0hat)
%
%  Atilde   =  Ctilde(uhat,t0hat)
%
%  Bhat     =  Chat(uhat,t0tilde)
%
%  Btilde   =  Ctilde(uhat,t0tilde)
%
%  Dtilde   =  Ctilde(t0tilde,t0tilde)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t3      =  t^3-term
%
%  Ehat    =  Chat(t0hat,t0hat)
%
%  Etilde  =  Ctilde(t0hat,t0hat)
%
%  Fhat    =  Chat(t0hat,t0tilde)
%
%  Ftilde  =  Ctilde(t0hat,t0tilde)

E = convolution_sum_Burgers(t0hat,t0hat,alpha);
Ehat = E;
Etilde = E;
Ehat(G_modes) = 0;
Etilde(F_modes) = 0;

F = convolution_sum_Burgers(t0hat,t0tilde,alpha);
Fhat = F;
Ftilde = F;
Fhat(G_modes) = 0;
Ftilde(F_modes) = 0;

int1 = convolution_sum_Burgers(u_full,Ahat - 2*Bhat - 2*Atilde + Btilde,alpha);
int1(F_modes) = 0;

term1 = convolution_sum_Burgers(u_full,2*int1 + Etilde - Ftilde + Dtilde,alpha);
term1(G_modes) = 0;

term2 = convolution_sum_Burgers(t0tilde,-Atilde + Btilde,alpha);
term2(G_modes) = 0;


t3 = 4*term1 + 12*term2;