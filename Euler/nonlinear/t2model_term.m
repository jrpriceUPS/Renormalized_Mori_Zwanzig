function [t2,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2)
%
% Computes the RHS for every mode in the full model for 3D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    u_full  =  full array of current Fourier state (2Mx2Mx2Mx3)
%
%     t0hat  =  full array of current Fourier state of C_hat(u,u)
%
%   t0tilde  =  full array of current Fourier state of C_tilde(u,u)
%
%     t1hat  =  full array of current Fourier state of hat{t1-term}
%
%         a  =  indices of positive resolved modes 1:M
%
%         b  =  indices of negative resolved modes -M:-1
%
%         k  =  array of wavenumbers (2Mx2Mx2Mx3)
%
%   a_tilde  =  indices of unresolved modes
%
%  a_tilde2  =  indices corresponding to modes included only for
%               dealiasing
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%       t2  =  t-model term of derivative of each resolved mode
%
%    t2hat  =  resolved part of the t-model term of derivative of each resolved mode
%
%  t2tilde  =  unresolved part of the t-model term of derivative of each resolved mode

[~,Ahat,Atilde] = Dk(u_full,t0hat,a,b,k,a_tilde,a_tilde2);
[B,Bhat,Btilde] = Ck(t0tilde,t0tilde,a,b,k,a_tilde,a_tilde2);

[first_term,~,~] = Dk(u_full,t1tilde-Atilde,a,b,k,a_tilde,a_tilde2);

t2 = 1/2*first_term + B;