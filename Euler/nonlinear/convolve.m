function convolution = convolve(u_full,v_full)
%
% Computes a convolution sum of two quantities by transforming them into
% real space, computing a product, and transforming back into Fourier
% space. The quantity it computes is:
%
% convolution_k = sum_{p+q = k} u_p v_q
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
% u_full  =  an array in the convolution (2M x 2M x 2M x 3)
%
% v_full  =  an array in the convolution (2M x 2M x 2M x 3)
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% convolution  =  the convolution sum of u and v (2M x 2M x 2M x 2M x 3)

u_real = ifftn_norm(u_full);
v_real = ifftn_norm(v_full);

s = size(u_real);

uv_T = zeros([s 3]);

for i = 1:s(1)
    for j = 1:s(2)
        for k = 1:s(3)
            
            uv_T(i,j,k,:,:) = squeeze(u_real(i,j,k,:))*squeeze(v_real(i,j,k,:)).';
            
        end
    end
end

convolution = fftn_norm_conv(uv_T);