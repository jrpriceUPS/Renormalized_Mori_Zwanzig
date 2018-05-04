function du_dt = markov_RHS(params)
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
%        N  =  maximum positive wavenumber
%
%  a_tilde  =  indices of unresolved modes
%
% a_tilde2  =  indices corresponding to modes included only for
%               dealiasing
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
a_tilde2 = params.a_tilde2;

% the full model is a simple convolution Ck(u,u)
du_dt = markov_term(u_full,a,b,k,a_tilde,a_tilde2);

du_dt = u_squishify(du_dt,N);