function [nonlin1,uu_star] = tmodel_term_BKdV(u_full,nonlin0,F_modes)
%
%Computes the t-model term of Burgers-KdV for a given state vector
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
%  alpha    =  coefficient on the nonlinear term
%
%  F_modes  =  vector of which modes in u_full correspond to resolved modes
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  nonlin1  =  the t-model term
%
%  uu_star  =  the unresolved part of the markov convolution

%eliminate resolved modes of first convolution
uu_star = nonlin0;
uu_star(F_modes)=0;

%compute t model term
nonlin1 = 2*convolution_sum_BKdV(u_full,uu_star);