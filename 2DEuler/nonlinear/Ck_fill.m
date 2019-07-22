function [C] = Ck_fill(irange,jrange,k,convo)
%
% A function to use the computed convolution of v and w to fill out the
% array for Ck(v,w)
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%  irange  =  the range of k_x values (positive, negative, or zero)
%
%  jrange  =  the range of k_y values (positive, negative, or zero)
%
%       k  =  the array of wavevectors (2Mx2Mx2)
%
%   convo  =  the convolution sum ( sum_{p+q=k} v_q*w_p^T )
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  C  =  the convolution C_k(v,w) for each Fourier mode in the range
%            [irange,jrange,:]

% find the size of the arrays we are using
s = size(convo);


t1 = zeros(s(1:end-1));
k1 = repmat(k(irange,jrange,1),[1,1,2]);
k2 = repmat(k(irange,jrange,2),[1,1,2]);
t1(irange,jrange,:) = -1j*(k1.*convo(irange,jrange,:,1) + k2.*convo(irange,jrange,:,2));


t2 = zeros(s(1:end-1));

wave_vec_mag = sum(k(irange,jrange,:).^2,3);
wave_vec_mag(wave_vec_mag==0)=1;
wave_vec_mag = repmat(wave_vec_mag,[1,1,2]);

t1_1 = repmat(t1(irange,jrange,1),[1,1,2]);
t1_2 = repmat(t1(irange,jrange,2),[1,1,2]);

B = k1.*t1_1 + k2.*t1_2;
t2(irange,jrange,:) = k(irange,jrange,:).*B./wave_vec_mag;



% retain only relevant terms
C = t1(irange,jrange,:)-t2(irange,jrange,:);
