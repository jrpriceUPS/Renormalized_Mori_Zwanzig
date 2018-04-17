function [nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon)
%
%Computes the BCH t^2-model term for a given state vector
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
%  alpha    =  coefficient on the nonlinear term
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%  G_modes  =  vector of which modes in u_full correspond to unresolved
%              modes
%
%  k        =  vector of wavenumbers corresponding to entries of u_full
%
%  epsilon  =  degree of dispersion
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin2       =  t^2-term
%
%  uk            =  u_full . * k (each entry is ku)
%
%  uu            =  resolved part of the markov convolution
%
%  uk_uu_u       =  resolved part of convolution of uu+uk with u_full
%
%  uk_uu_u_star  =  unresolved part of convolution of uu+uk with u_full



%eliminate unresolved modes of first convolution
uu = nonlin0;
uu(G_modes)=0;

%compute u_full with dispersive term
uk = 1j*epsilon^2*k.^3.*u_full;

%compute t^2 model term
uk_uu_u = convolution_sum(uk + uu,u_full,alpha);
uk_uu_u_star = uk_uu_u;
uk_uu_u(G_modes) = 0;
uk_uu_u_star(F_modes) = 0;

nonlin2 =  2*convolution_sum(uk+uu,uu_star,alpha) ...
         + 4*convolution_sum(u_full,uk_uu_u_star,alpha);