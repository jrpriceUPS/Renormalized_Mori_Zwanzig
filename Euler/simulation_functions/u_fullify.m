function u_full = u_fullify(u,M)
%
% Takes a compressed state array and uses the appropriate complex
% conjugates to construct the full state array
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
% u  =  an NxNxNx3x4 array containing the Fourier coefficients for
%       half the necessary modes
%
% M  =  the desired size of the full system
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% u_full  =  a 2Mx2Mx2Mx3 array containing the Fourier coefficients for all
%            Fourier modes


s = size(u); % find size of array
N = s(1);

% create hotkeys for different index ranges
a = 2:N;
b = 2*M:-1:2*M-N+2;

u_full = zeros(2*M,2*M,2*M,3); % create output array

% fill in values from reduced array
% we save data for wavevectors k = [ , , ]
% let + mean positive entries in that coordinate of k,
%     - means negative entries in that coordinate of k
%     0 means zeros in that coordinate of k


u_full(a,a,a,:) = u(a,a,a,:,1); % k = [+,+,+]
u_full(a,a,b,:) = u(a,a,a,:,2); % k = [+,+,-]
u_full(a,b,a,:) = u(a,a,a,:,3); % k = [+,-,+]
u_full(b,a,a,:) = u(a,a,a,:,4); % k = [-,+,+]

u_full(1,a,a,:) = u(1,a,a,:,1); % k = [0,+,+]
u_full(a,1,a,:) = u(a,1,a,:,1); % k = [+,0,+]
u_full(a,a,1,:) = u(a,a,1,:,1); % k = [+,+,0]

u_full(1,a,b,:) = u(1,a,a,:,2); % k = [0,+,-]
u_full(a,1,b,:) = u(a,1,a,:,2); % k = [+,0,-]
u_full(a,b,1,:) = u(a,a,1,:,2); % k = [+,-,0]

u_full(1,1,a,:) = u(1,1,a,:,1); % k = [0,0,+]
u_full(1,a,1,:) = u(1,a,1,:,1); % k = [0,+,0]
u_full(a,1,1,:) = u(a,1,1,:,1); % k = [+,0,0]

u_full(1,1,1,:) = u(1,1,1,:,1); % k = [0,0,0]



% fill in related conjugate values
u_full(b,b,b,:) = conj(u_full(a,a,a,:)); % k = [-,-,-]
u_full(b,b,a,:) = conj(u_full(a,a,b,:)); % k = [-,-,+]
u_full(b,a,b,:) = conj(u_full(a,b,a,:)); % k = [-,+,-]
u_full(a,b,b,:) = conj(u_full(b,a,a,:)); % k = [+,-,-]

u_full(1,b,b,:) = conj(u_full(1,a,a,:)); % k = [0,-,-]
u_full(b,1,b,:) = conj(u_full(a,1,a,:)); % k = [-,0,-]
u_full(b,b,1,:) = conj(u_full(a,a,1,:)); % k = [-,-,0]

u_full(1,b,a,:) = conj(u_full(1,a,b,:)); % k = [0,-,+]
u_full(b,1,a,:) = conj(u_full(a,1,b,:)); % k = [-,0,+]
u_full(b,a,1,:) = conj(u_full(a,b,1,:)); % k = [-,+,0]

u_full(b,1,1,:) = conj(u_full(a,1,1,:)); % k = [-,0,0]
u_full(1,b,1,:) = conj(u_full(1,a,1,:)); % k = [0,-,0]
u_full(1,1,b,:) = conj(u_full(1,1,a,:)); % k = [0,0,-]