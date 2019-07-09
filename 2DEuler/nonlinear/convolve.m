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
% u_full  =  an array in the convolution (2M x 2M x 2)
%
% v_full  =  an array in the convolution (2M x 2M x 2)
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% convolution  =  the convolution sum of u and v (2M x 2M x 2M x 2)
% 
%

u_real = ifftn_norm(u_full);
v_real = ifftn_norm(v_full);

s = size(u_real);

v_1 = repmat(v_real(:,:,1),[1,1,2]);
v_2 = repmat(v_real(:,:,2),[1,1,2]);

uv_T = zeros([s 2]);
uv_T(:,:,:,1) = u_real.*v_1;
uv_T(:,:,:,2) = u_real.*v_2;


convolution = fftn_norm_conv(uv_T);
