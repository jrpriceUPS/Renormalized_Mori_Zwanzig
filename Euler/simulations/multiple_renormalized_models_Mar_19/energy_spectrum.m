function [spectrum,k_list] = energy_spectrum(u_full,M)

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

k_mag = sum(k.^2,4);

u_mag = sqrt(sum(u_full.*conj(u_full),4));

k_list = 1:M/3;
spectrum = zeros(length(k_list),1);

for i = 1:length(k_list)
    current_k = k_list(i);
    spectrum(i) = sum(u_mag(current_k-1/2 < k_mag & k_mag < current_k+1/2));
end