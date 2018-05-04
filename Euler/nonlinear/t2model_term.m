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
%   t1tilde  =  full array of current Fourier state of tilde{t1-term}
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
%       t2  =  t^2-model term of derivative of each resolved mode
%
%     Ahat  =  resolved part of Dk(hat{u},hat{t0})
%
%   Atilde  =  unresolved part of Dk(hat{u},hat{t0})
%
%     Bhat  =  resolved part of Dk(hat{u},tilde{t1}-Atilde)
%
%   Btilde  =  unresolved part of Dk(hat{u},tilde{t1}-Atilde)

[~,Ahat,Atilde] = Dk(u_full,t0hat,a,b,k,a_tilde,a_tilde2);
[B,Bhat,Btilde] = Ck(t0tilde,t0tilde,a,b,k,a_tilde,a_tilde2);

[first_term,~,~] = Dk(u_full,t1tilde-Atilde,a,b,k,a_tilde,a_tilde2);

t2 = -first_term - 2*B;