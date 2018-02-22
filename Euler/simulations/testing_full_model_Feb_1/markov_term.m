function [du_dt,du_dt_hat,du_dt_tilde] = markov_term(u_full,a,b,k,a_hat,a_tilde)
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
%        a  =  indices of positive resolved modes 1:M
%
%        b  =  indices of negative resolved modes -M:-1
%
%        k  =  array of wavenumbers (2Mx2Mx2Mx3)
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  du_dt  =  hello

% the full model is a simple convolution Ck(u,u)
[du_dt,du_dt_hat,du_dt_tilde] = Ck(u_full,u_full,a,b,k,a_hat,a_tilde);