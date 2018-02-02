function helicity = helicity(u_array)
%
% Computes the total helicity at all times for an output array. This is
% defined as:
%
%  w  =  int u dot (grad cross u) dx
%
%  We calculate grad cross u in fourier space, and the dot product in real
%  space
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
%  helicity  =  a length(t)x1 array of the helicity at all times
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


helicity = zeros(s(end),1);

% loop through all times
for i = 1:s(end)
    
    % isolate the current state and create permutations for curl
    % computation
    u = u_array(:,:,:,:,:,i);
    u_full1 = u_fullify(u,N);
    u_full2 = u_full1(:,:,:,[2,3,1]);
    u_full3 = u_full1(:,:,:,[3,1,2]);
    
    % calculate the curl in Fourier space and transform it into real space
    curl_fourier = 1j*(k2.*u_full3 - k3.*u_full2);
    curl = ifftn_norm(curl_fourier);
    
    % calculate the real space u and then compute the helicity for the
    % given time
    u_real = ifftn_norm(u_full1);
    helicity(i) = sum(sum(sum(sum(u_real.*curl))));
    
end