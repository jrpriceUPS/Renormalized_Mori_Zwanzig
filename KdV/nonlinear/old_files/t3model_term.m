function nonlin3 = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star)
%
%Computes the BCH t^3-model term for a given state vector
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
%  k        =  vector of wavenumbers corresponding to entries of u_full
%
%  epsilon  =  degree of dispersion
%
%  u_full   =  a full state vector (positive and negative modes)
%
%  uu       =  resolved part of the markov convolution
%
%  uu_star  =  the resolved part of the markov convolution
%
%  uk       =  the resolved part of k^3 .* u_full
%
%  uk_uu_u       =  resolved part of convolution of uu+uk with u_full
%
%  uk_uu_u_star  =  unresolved part of convolution of uu+uk with u_full
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin3  =  t^3-term

%compute auxiliary convolutions
k_uk_uu = 1j*epsilon^2*k.^3.*(uk+uu);

third_term = convolution_sum(k_uk_uu+2*uk_uu_u,u_full,alpha);
third_term(F_modes) = 0;

fourth_term = convolution_sum(uk+uu,uk+uu,alpha);
fourth_term(F_modes) = 0;


%compute BCH 3rd order term
nonlin3 =  2*convolution_sum(k_uk_uu+2*uk_uu_u,uu_star,alpha) ...
         + 8*convolution_sum(uk+uu,uk_uu_u_star,alpha) ...
         + 4*convolution_sum(u_full,third_term,alpha) ...
         + 4*convolution_sum(u_full,fourth_term,alpha);
     
