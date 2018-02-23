function [t3,t3hat,t3tilde,E,Ehat,Etilde,F,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde)
%
% Computes the RHS for every mode in the t^3-model term for 3D Euler
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%   u_full  =  full array of current Fourier state (2Mx2Mx2Mx3)
%
%    t0hat  =  full array of current Fourier state of C_hat(u,u)
%
%  t0tilde  =  full array of current Fourier state of C_tilde(u,u)
%
%    t1hat  =  full array of current Fourier state of hat{t1-term}
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
%       t3  =  t^3-model term of derivative of each resolved mode
%
%    t3hat  =  resolved part of the t^3-model term of derivative of each resolved mode
%
%  t3tilde  =  unresolved part of the t^3-model term of derivative of each resolved mode


[E,Ehat,Etilde] = Dk(t0tilde,t0hat,a,b,k,a_tilde);
[F,Fhat,Ftilde] = Ck(t0hat,t0hat,a,b,k,a_tilde);

[term1,t1_hat,t1_tilde] = Dk(t0tilde,t1tilde-Atilde,a,b,k,a_tilde);
[~,~,t2_tilde] = Dk(u_full,Ahat-2*t1hat+t1tilde-2*Atilde,a,b,k,a_tilde);
[term3,t3_hat,t3_tilde] = Dk(u_full,t2_tilde+2*Btilde-Etilde+2*Ftilde,a,b,k,a_tilde);


t3 = 2*term1 + 1/6*term3;
t3hat = 2*t1_hat + 1/6*t3_hat;
t3tilde = 2*t1_tilde + 1/6*t3_tilde;