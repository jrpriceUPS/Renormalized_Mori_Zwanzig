function [C,Chat,Ctilde] = Ck(v_full,w_full,a,b,k,a_tilde,a_tilde2)
%
% Computes the convolution of v and w as well as the resolved and
% unresolved versions of the same
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
%        C  =  Ck(v,w) array in compressed form (MxMx2x2)
%
%    C_hat  =  C_hat(v,w) array in compressed form (MxMx2x2)
%
%  C_tilde  =  C_tilde(v,w) array in compressed form (MxMx2x2)


% begin by computing the challenging sum part as an outer product in real
% space
convo = convolve(v_full,w_full);

% identify the size of the array and construct the output
M = a(end);
C = zeros(M,M,2,2);


% fill in the modes in a compressed form, keeping only the ones that matter
C(a,a,:,1) = Ck_fill(a,a,k,convo);
C(a,a,:,2) = Ck_fill(a,b,k,convo);

C(1,a,:,1) = Ck_fill(1,a,k,convo);
C(a,1,:,1) = Ck_fill(a,1,k,convo);

C(1,a,:,2) = Ck_fill(1,b,k,convo);

C(1,1,:,1) = Ck_fill(1,1,k,convo);


% clear out modes for resolved and unresolved parts of the array
Chat = mode_clearer(C,a_tilde);
Ctilde = C - Chat;
Ctilde = mode_clearer(Ctilde,a_tilde2);

C = u_fullify(C,M);
Chat = u_fullify(Chat,M);
Ctilde = u_fullify(Ctilde,M);

