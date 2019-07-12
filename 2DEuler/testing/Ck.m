function [C1,C2] = Ck(v_full,w_full,a,b,k,a_tilde,a_tilde2)
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
C1 = zeros(M,M,2,2);
C2 = zeros(M,M,2,2);

% fill in the modes in a compressed form, keeping only the ones that matter
[C1(a,a,:,1),C2(a,a,:,1)] = Ck_fill(a,a,k,convo);
[C1(a,a,:,2),C2(a,a,:,2)] = Ck_fill(a,b,k,convo);

[C1(1,a,:,1),C2(1,a,:,1)] = Ck_fill(1,a,k,convo);
[C1(a,1,:,1),C2(a,1,:,1)] = Ck_fill(a,1,k,convo);

[C1(1,a,:,2),C2(1,a,:,2)] = Ck_fill(1,b,k,convo);

[C1(1,1,:,1),C2(1,1,:,1)] = Ck_fill(1,1,k,convo);


% clear out modes for resolved and unresolved parts of the array
C1hat = mode_clearer(C1,a_tilde);
C1tilde = C1 - C1hat;
C1tilde = mode_clearer(C1tilde,a_tilde2);

C1 = u_fullify(C1,M);
C1hat = u_fullify(C1hat,M);
C1tilde = u_fullify(C1tilde,M);

C2hat = mode_clearer(C2,a_tilde);
C2tilde = C2 - C2hat;
C2tilde = mode_clearer(C2tilde,a_tilde2);

C2 = u_fullify(C2,M);
C2hat = u_fullify(C2hat,M);
C2tilde = u_fullify(C2tilde,M);

