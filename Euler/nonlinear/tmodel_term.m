function [t1,t1hat,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde)
%
% Computes the RHS for every mode in the full model for 3D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%   u_full  =  full array of current Fourier state (2Mx2Mx2Mx3)
%
%  t0tilde  =  full array of current Fourier state of C_tilde(u,u)
%
%        a  =  indices of positive resolved modes 1:M
%
%        b  =  indices of negative resolved modes -M:-1
%
%        k  =  array of wavenumbers (2Mx2Mx2Mx3)
%
%  a_tilde  =  indices of unresolved modes
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  t1  =  t-model term of derivative of each resolved mode

% the t-model is Dk(u_hat,C_tilde(u_hat,u_hat))
[t1,t1hat,t1tilde] = Dk(u_full,t0tilde,a,b,k,a_tilde);