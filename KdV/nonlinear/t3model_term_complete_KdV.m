function [nonlin3,uk6,E,E_star,F,F_star] = t3model_term_complete_KdV(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk3,A,A_star,B,B_star,C,C_star,D_star)
%
%Computes the complete t^3-model term of KdV for a given state vector
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
%  k        =  vector of wavenumbers corresponding to entries of u_full
%
%  epsilon  =  degree of dispersion
%
%  u_full   =  a full state vector (positive and negative modes)
%
%  uu       =  resolved part of the markov convolution (entries: C_k(u,u))
%
%  uu_star  =  unresolved part of the markov convolution (entries: C^*_k(u,u))
%
%  uk3      =  u_full . * k.^3 (entries: k^3*u_k)
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
%  D_star   =  unresolved part of the convolution of uu_star with itself
%              (entries: C^*_k(C^*(u,u),C^*(u,u)))
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin3  =  t^3-term
%
%  uk6      =  u_full . * k.^6 (entries: k^6*u_k)
%
%  E        =  resolved part of convolution of i*epsilon^2*uk3+uu with
%              itself (entries: C_k(i*epsilon^2*k^3*u+C(u,u),i*epsilon^2*k^3*u+C(u,u)))
%
%  E_star   =  unresolved part of convolution of i*epsilon^2*uk3+uu with
%              itself (entries: C^*_k(i*epsilon^2*k^3*u+C(u,u),i*epsilon^2*k^3*u+C(u,u)))
%
%  F        =  resolved part of convolution of i*epsilon^2*uk3+uu with uu_star
%              (entries: C_k(i*epsilon^2*k^3*u+C(u,u),C^*(u,u)))
%
%  F_star   =  unresolved part of convolution of i*epsilon^2*uk3+uu with uu_star
%              (entries: C^*_k(i*epsilon^2*k^3*u+C(u,u),C^*(u,u)))

%compute complete t^3-model term (and auxiliary terms)
uk6 = k.^3.*uk3;

E = convolution_sum_KdV(1j*epsilon^2*uk3+uu,1j*epsilon^2*uk3+uu,alpha);
E_star = E;
E(G_modes) = 0;
E_star(F_modes) = 0;

F = convolution_sum_KdV(uu_star,1j*epsilon^2*uk3+uu,alpha);
F_star = F;
F(G_modes) = 0;
F_star(F_modes) = 0;


int1 = -2*B_star+C_star;
int2 = convolution_sum_KdV(u_full,-epsilon^4.*uk6+1j*epsilon^2.*(A+A_star)...
                              +2*(B-2*C)+2*(C_star-2*B_star),alpha);
int2(F_modes) = 0;
int3 = E_star-F_star;
int4 = D_star;
int5 = C_star - B_star;

nonlin3 = 2*convolution_sum_KdV(u_full,-k.^3*epsilon^4.*A_star + 2*1j*epsilon^2*k.^3.*int1...
                            + 2*int2 + 2*int3 + 2*int4,alpha)...
          +6*convolution_sum_KdV(uu_star,1j*epsilon^2.*A_star + 2*int5,alpha);