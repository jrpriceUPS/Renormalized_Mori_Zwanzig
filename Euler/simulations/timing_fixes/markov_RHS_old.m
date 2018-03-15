function du_dt = markov_RHS_old(params)
%
% Computes the RHS for every mode in the markov model for 3D Euler
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
%  a_tilde  =  indices of unresolved modes
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  du_dt  =  derivative of each mode

u_full = params.u_full;
a = params.a;
b = params.b;
k = params.k;
N = params.N;
a_tilde = params.a_tilde;

% the full model is a simple convolution Ck(u,u)
du_dt = markov_term_old(u_full,a,b,k,a_tilde);

du_dt = u_squishify(du_dt,N);