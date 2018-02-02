function normalized_ifftn = ifftn_norm(u_full)
%
%Computes the n-dimensional IFFT of the NxNxN u_full using the more standard
%normalization of 1/N^3.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_full  =  Fourier space vector to be transformed
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  normalized_ifftn  =  normalized IFFT of that vector


normalized_ifftn = zeros(size(u_full));
for i = 1:3
    normalized_ifftn(:,:,:,i) = real(ifftn(u_full(:,:,:,i))*numel(u_full(:,:,:,i)));
end