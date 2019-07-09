function [D,Dhat,Dtilde] = Dk(v_full,w_full,a,b,k,a_tilde,a_tilde2)
%
% Computes the double convolution Dk(v,w) = Ck(v,w) + Ck(w,v) and the
% associated resolved Dhat(vw) and unresolved Dtilde(v,w) versions
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%    v_full  =  full Fourier form of first argument of C (2Mx2Mx2)
%
%    w_full  =  full Fourier form of second argument of C (2Mx2Mx2)
%
%         a  =  indices corresponding to positive modes 1:M
%
%         b  =  indices corresponding to negative modes -M:-1
%
%         k  =  array of wavevectors (2Mx2Mx2)
%
%   a_tilde  =  indices corresponding to positive unresolved
%
%  a_tilde2  =  indices corresponding to modes included only for
%               dealiasing
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%        D  =  Dk(v,w) array in compressed form (MxMx2x2)
%
%    D_hat  =  D_hat(v,w) array in compressed form (MxMx2x2)
%
%  D_tilde  =  D_tilde(v,w) array in compressed form (MxMx2x2)

[C1,C1hat,C1tilde] = Ck(v_full,w_full,a,b,k,a_tilde,a_tilde2);
[C2,C2hat,C2tilde] = Ck(w_full,v_full,a,b,k,a_tilde,a_tilde2);

D = C1 + C2;
Dhat = C1hat + C2hat;
Dtilde = C1tilde + C2tilde;