function nonlin4 = t4model_term_complete(alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk3,uk6,A,A_star,B,B_star,C,C_star,D,D_star,E,E_star,F,F_star)
%
%Computes the complete t^4-model term for a given state vector
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
%  uk6      =  u_full . * k.^6 (entries: k^6*u_k)
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
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin4  =  t^4-term

%compute complete t^4-model term (and auxiliary terms)

internal1 = convolution_sum(u_full,-epsilon^4*uk6+1j*epsilon^2*(A+A_star)...
                                   +2*B - 4*C - 4*B_star + 2*C_star,alpha);
internal1(F_modes) = 0;

internal2 = 1j*epsilon^2*k.^3.*convolution_sum(u_full,-3*epsilon^4*uk6 ...
                                                      +1j*epsilon^2*(3*A+A_star) ...
                                                      -2*(-3*B+5*C) ...
                                                      +2*(-3*B_star+C_star),alpha);
internal2(F_modes) = 0;                      
                          
auxiliary1 = 2*convolution_sum(u_full,epsilon^4*uk6-1j*epsilon^2*(A+3*A_star)...
                                      +2*(3*C-B)+2*(5*B_star-3*C_star),alpha);
auxiliary1(G_modes) = 0;
auxiliary2 = 2*convolution_sum(u_full,-3*epsilon^4*uk6+1j*epsilon^2*(3*A+A_star)...
                                      +2*(3*B-5*C)+2*(-3*B_star+C_star),alpha);
auxiliary2(F_modes) = 0;
internal3 = convolution_sum(u_full,1j*k.^3.*uk6*epsilon^6 ...
                                   +k.^3*epsilon^4.*(A-A_star) ...
                                   +2*1j*epsilon^2*k.^3.*(3*C-B) ...
                                   +2*1j*epsilon^2*k.^3.*(-3*B_star+C_star) ...
                                   +auxiliary1 ...
                                   +auxiliary2 ...
                                   -2*(E-2*F) ...
                                   +2*(3*E_star-2*F_star) ...
                                   -6*D + 2*D_star,alpha);
internal3(F_modes) = 0;
                          
internal4 = convolution_sum(1j*epsilon^2*uk3+uu,3*epsilon^4*uk6 - 1j*epsilon^2*(3*A+A_star) ...
                                                  +2*(-3*B+5*C)+2*(3*B_star-C_star),alpha);
internal4(F_modes) = 0;

internal5 = convolution_sum(uu_star,-epsilon^4*uk6 ...
                                    +1j*epsilon^2*(A+3*A_star) ...
                                    +2*B - 6*C - 10*B_star + 6*C_star,alpha);
internal5(F_modes) = 0;

nonlin4 = 2*convolution_sum(u_full,-1j*epsilon^6*k.^6.*A_star...
                                   +2*k.^6*epsilon^4.*(3*B_star-C_star)...
                                   +2*internal2 ...
                                   +2*internal3 ...
                                   +2*internal4 ...
                                   -2*k.^3*1j*epsilon^2.*(2*F_star-3*E_star) ...
                                   +2*k.^3*1j*epsilon^2.*D_star ...
                                   +2*internal5,alpha) ...
         +8*convolution_sum(uu_star,-k.^3*epsilon^4.*A_star ...
                                    +2*1j*epsilon^2*k.^3.*(-2*B_star+C_star)...
                                    +2*internal1...
                                    +2*(E_star-F_star)...
                                    +2*D_star,alpha)...
         -48*convolution_sum(B_star,1j*epsilon^2*A_star + 2*C_star,alpha)...
         +6*convolution_sum(1j*epsilon^2*A_star+2*(B_star+C_star),...
                            1j*epsilon^2*A_star+2*(B_star+C_star),alpha);
                        

nonlin4 = -nonlin4;