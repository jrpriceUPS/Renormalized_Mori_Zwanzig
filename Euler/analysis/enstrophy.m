function enstrophy = enstrophy(u_array)
%
% Computes the enstrophy at all times for an output array. This is
% defined as:
%
%  w  =  int |grad cross u|^2 dx
%
%  We calculate grad cross u in fourier space, and use parseval's equality
%  to compute the integral from those coefficients
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%  u_array  =  an NxNxNx3x4xlength(t) array of Fourier states
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  enstrophy  =  a length(t)x1 array of the enstrophy at all times


s = size(u_array);
N = s(1);

% make k array
k_vec = [0:N-1,-N:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k2 = zeros(2*N,2*N,2*N,3);
k3 = zeros(2*N,2*N,2*N,3);

% create permutations for curl computation
k2(:,:,:,1) = ky;
k2(:,:,:,2) = kz;
k2(:,:,:,3) = kx;

k3(:,:,:,1) = kz;
k3(:,:,:,2) = kx;
k3(:,:,:,3) = ky;


enstrophy = zeros(s(end),1);

% loop through all times
for i = 1:s(end)
    
    % isolate the current state and create permutations for curl
    % computation
    u = u_array(:,:,:,:,:,i);
    u_full1 = u_fullify(u,N);
    u_full2 = u_full1(:,:,:,[2,3,1]);
    u_full3 = u_full1(:,:,:,[3,1,2]);
    
    % calculate the curl in Fourier space
    curl = 1j*(k2.*u_full3 - k3.*u_full2);
    enstrophy(i) = sum(curl(:).*conj(curl(:)));
    
end