function u = u_squishify(u_full,N)
%
%Reduces a full array in Fourier space to a condensed version that does not
%duplicate conjugate data
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%  u_full  =  2Mx2Mx2 fully constructed Fourier coefficient array
%
%       N  =  maximum value of wavevectors to include in squished version
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  u  =  NxNx2x2 array of only the non-conjugate Fourier modes (cuts
%        array size in half or more)

% compute number of modes in the original array
s = size(u_full);
M = s(1)/2;

% construct output
u = zeros(N,N,2,2);

% a represents positive modes, while b is negative modes (listed in
% reverse)
a = 2:N;
b = 2*M:-1:2*M-N+2;

% fill in the modes we intend to keep
% in each case, we save data for wavevectors k = [ , ]
% let + mean positive entries in that coordinate of k,
%     - means negative entries in that coordinate of k
%     0 means zeros in that coordinate of k

u(a,a,:,1) = u_full(a,a,:); % k = [+,+] (conjugate is [-,-])
u(a,a,:,2) = u_full(a,b,:); % k = [+,-] (conjugate is [-,+])

u(1,a,:,1) = u_full(1,a,:); % k = [0,+] (conjugate is [0,-])
u(a,1,:,1) = u_full(a,1,:); % k = [+,0] (conjugate is [-,0])

u(1,1,:,1) = u_full(1,1,:); % k = [0,0]