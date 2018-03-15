function [C,Chat,Ctilde] = Ck_old(v_full,w_full,a,b,k,a_tilde)
%
% Computes the convolution of v and w as well as the resolved and
% unresolved versions of the same
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%   v_full  =  full Fourier form of first argument of C (2Mx2Mx2Mx3)
%
%   w_full  =  full Fourier form of second argument of C (2Mx2Mx2Mx3)
%
%        a  =  indices corresponding to positive modes 1:M
%
%        b  =  indices corresponding to negative modes -M:-1
%
%        k  =  array of wavevectors (2Mx2Mx2Mx3)
%
%  a_tilde  =  indices corresponding to positive unresolved
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%        C  =  Ck(v,w) array in compressed form (MxMxMx3x4)
%
%    C_hat  =  C_hat(v,w) array in compressed form (MxMxMx3x4)
%
%  C_tilde  =  C_tilde(v,w) array in compressed form (MxMxMx3x4)


% begin by computing the challenging sum part as an outer product in real
% space
convo = convolve_old(v_full,w_full);

% identify the size of the array and construct the output
M = a(end);
C = zeros(M,M,M,3,4);


% fill in the modes in a compressed form, keeping only the ones that matter
C(a,a,a,:,1) = Ck_fill_old(a,a,a,k,convo);
C(a,a,a,:,2) = Ck_fill_old(a,a,b,k,convo);
C(a,a,a,:,3) = Ck_fill_old(a,b,a,k,convo);
C(a,a,a,:,4) = Ck_fill_old(b,a,a,k,convo);

C(1,a,a,:,1) = Ck_fill_old(1,a,a,k,convo);
C(a,1,a,:,1) = Ck_fill_old(a,1,a,k,convo);
C(a,a,1,:,1) = Ck_fill_old(a,a,1,k,convo);

C(1,a,a,:,2) = Ck_fill_old(1,a,b,k,convo);
C(a,1,a,:,2) = Ck_fill_old(a,1,b,k,convo);

C(1,a,a,:,3) = Ck_fill_old(1,b,a,k,convo);
C(a,a,1,:,3) = Ck_fill_old(a,b,1,k,convo);

C(a,1,a,:,4) = Ck_fill_old(b,1,a,k,convo);
C(a,a,1,:,4) = Ck_fill_old(b,a,1,k,convo);

C(1,1,a,:,1) = Ck_fill_old(1,1,a,k,convo);
C(1,a,1,:,1) = Ck_fill_old(1,a,1,k,convo);
C(a,1,1,:,1) = Ck_fill_old(a,1,1,k,convo);

C(1,1,1,:,1) = Ck_fill_old(1,1,1,k,convo);


% clear out modes for resolved and unresolved parts of the array
Chat = mode_clearer(C,a_tilde);
Ctilde = C - Chat;

C = u_fullify(C,M);
Chat = u_fullify(Chat,M);
Ctilde = u_fullify(Ctilde,M);

