function u_full = u_fullify(u,M)
%
% Takes a compressed state array and uses the appropriate complex
% conjugates to construct the full state array
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
% u  =  an NxNx2x2 array containing the Fourier coefficients for
%       half the necessary modes
%
% M  =  the desired size of the full system
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% u_full  =  a 2Mx2Mx2 array containing the Fourier coefficients for all
%            Fourier modes


s = size(u); % find size of array
N = s(1);

% create hotkeys for different index ranges
a = 2:N;
b = 2*M:-1:2*M-N+2;

u_full = zeros(2*M,2*M,2); % create output array

% fill in values from reduced array
% we save data for wavevectors k = [ , ]
% let + mean positive entries in that coordinate of k,
%     - means negative entries in that coordinate of k
%     0 means zeros in that coordinate of k

u_full(a,a,:) = u(a,a,:,1); % k = [+,+]
u_full(a,b,:) = u(a,a,:,2); % k = [+,-]

u_full(1,a,:) = u(1,a,:,1); % k = [0,+]
u_full(a,1,:) = u(a,1,:,1); % k = [+,0]

u_full(1,1,:) = u(1,1,:,1); % k = [0,0]



% fill in related conjugate values
u_full(b,b,:) = conj(u_full(a,a,:)); % k = [-,-]
u_full(b,a,:) = conj(u_full(a,b,:)); % k = [-,+]

u_full(b,1,:) = conj(u_full(a,1,:)); % k = [-,0]
u_full(1,b,:) = conj(u_full(1,a,:)); % k = [0,-]


