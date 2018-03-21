function [vort1,vort2] = vorticity(u_array)
%
% Computes the maximum of the vorticity (in infinity norm and two norm)
% at all times for an output array. This is defined as:
%
%  vort1  =  max(curl(u),inf)
%  vort2  =  max(curl(u),2)
%
%  We calculate grad cross u in fourier space, then transform back to real
%  space and compute both maxima
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


vort1 = zeros(s(end),1);
vort2 = zeros(s(end),1);

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
    real_curl = ifftn_norm(curl);
    vort1(i) = max(abs(real_curl(:)));
    vort_2norm = sqrt(sum(real_curl.^2,4));
    vort2(i) = max(vort_2norm(:));
    
end