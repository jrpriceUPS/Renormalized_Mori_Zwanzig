function du_dt = full_RHS(params)
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
%        N  =  maximum positive wavenumber
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  du_dt  =  derivative using the full model

u_full = params.u_full;
a = params.a;
b = params.b;
k = params.k;
N = params.N;

% the full model is a simple convolution Ck(u,u)
du_dt = Ck(u_full,u_full,a,b,k,[]);
du_dt = u_squishify(du_dt,N);