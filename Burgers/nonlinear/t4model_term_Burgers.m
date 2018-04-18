function t4 = t4model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde,Ehat,Etilde,Fhat,Ftilde)
%
%Computes the complete t^4-model term of KdV for a given state vector
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
%  Ehat    =  Chat(t0hat,t0hat)
%
%  Etilde  =  Ctilde(t0hat,t0hat)
%
%  Fhat    =  Chat(t0hat,t0tilde)
%
%  Ftilde  =  Ctilde(t0hat,t0tilde)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t4  =  t^4-term

%compute complete t^4-model term (and auxiliary terms)

int1a = convolution_sum_Burgers(u_full,2*(Ahat-3*Bhat)-10*Atilde+6*Btilde,alpha);
int1a(G_modes) = 0;

int1b = convolution_sum_Burgers(u_full,-6*Ahat+10*Bhat+2*(3*Atilde-Btilde),alpha);
int1b(F_modes) = 0;

int1 = convolution_sum_Burgers(u_full,int1a + Ehat - 2*Fhat + 3*Dhat + int1b - 3*Etilde + 2*Ftilde - Dtilde,alpha);
int1(F_modes) = 0;

int2 = convolution_sum_Burgers(t0hat,3*Ahat -5*Bhat - 3*Atilde + Btilde,alpha);
int2(F_modes) = 0;

int3 = convolution_sum_Burgers(t0tilde,-Ahat + 3*Bhat + 5*Atilde - 3*Btilde,alpha);
int3(F_modes) = 0;

term1 = convolution_sum_Burgers(u_full,int1 + int2 + int3,alpha);
term1(G_modes) = 0;


int4 = convolution_sum_Burgers(u_full,2*(-Ahat + 2*Bhat)+2*(2*Atilde - Btilde),alpha);
int4(F_modes) = 0;

term2 = convolution_sum_Burgers(t0tilde,int4 - Etilde + Ftilde - Dtilde,alpha);
term2(G_modes) = 0;


term3 = convolution_sum_Burgers(Atilde,Btilde,alpha);
term3(G_modes) = 0;


term4 = convolution_sum_Burgers(Atilde+Btilde,Atilde+Btilde,alpha);
term4(G_modes) = 0;


t4 = 8*term1 + 16*term2 + 96*term3 - 24*term4;
