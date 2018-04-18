function [nonlin2,uk3,uu,A,A_star,B,B_star,C,C_star,D,D_star] = t2model_term_complete_KdV(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon)
%
%Computes the complete t^2-model term of KdV for a given state vector
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
%  nonlin2  =  t^2-term
%
%  uk3      =  u_full . * k.^3 (entries: k^3*u_k)
%
%  uu       =  resolved part of the markov convolution (entries: C_k(u,u))
%
%  A        =  resolved part of k.^3 .* uu (entries: k^3*C_k(u,u))
%
%  A_star   =  unresolved part of k.^3 .* uu (entries: k^3*C^*_k(u,u))
%
%  B        =  resolved part of convolution of i*epsilon^2*uk3+uu with u_full
%              (entries: C_k(i*epsilon^2*k^3*u+C(u,u),u))
%
%  B_star   =  unresolved part of convolution of i*epsilon^2*uk3+uu with u_full
%              (entries: C^*_k(i*epsilon^2*k^3*u+C(u,u),u))
%
%  C        =  resolved part of convolution of uu_star with u_full 
%              (entries: C_k(C^*(u,u),u))
%
%  C_star   =  unresolved part of convolution of uu_star with u_full 
%              (entries: C^*_k(C^*(u,u),u))
%
%  D        =  resolved part of the convolution of uu_star with itself
%              (entries: C_k(C^*(u,u),C^*(u,u)))
%
%  D_star   =  unresolved part of the convolution of uu_star with itself
%              (entries: C^*_k(C^*(u,u),C^*(u,u)))


%eliminate unresolved modes of first convolution
uu = nonlin0;
uu(G_modes)=0;

%compute u_full with dispersive term
uk3 = k.^3.*u_full;

%compute complete t^2 model term (and auxiliary terms)
A = k.^3.*uu;
A_star = k.^3.*uu_star;

B = convolution_sum_KdV(1j*epsilon^2*uk3+uu,u_full,alpha);
B_star = B;
B(G_modes) = 0;
B_star(F_modes) = 0;

C = convolution_sum_KdV(uu_star,u_full,alpha);
C_star = C;
C(G_modes) = 0;
C_star(F_modes) = 0;

D = convolution_sum_KdV(uu_star,uu_star,alpha);
D_star = D;
D(G_modes) = 0;
D_star(F_modes) = 0;

nonlin2 =  -2*convolution_sum_KdV(u_full,1j*epsilon^2*A_star - 2*B_star + 2*C_star,alpha) - 2*D;