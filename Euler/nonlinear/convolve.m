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

for i = 1:2M
    for j = 1:2M
        uv_T(i,j,:,:) = squeeze(u(i,j,:))*squeeze(v(i,j,:))'
    end
end

s = size(u_real);

v_1 = repmat(v_real(:,:,:,1),[1,1,1,3]);
v_2 = repmat(v_real(:,:,:,2),[1,1,1,3]);
v_3 = repmat(v_real(:,:,:,3),[1,1,1,3]);

uv_T = zeros([s 3]);
uv_T(:,:,:,:,1) = u_real.*v_1;
uv_T(:,:,:,:,2) = u_real.*v_2;
uv_T(:,:,:,:,3) = u_real.*v_3;


convolution = fftn_norm_conv(uv_T);
