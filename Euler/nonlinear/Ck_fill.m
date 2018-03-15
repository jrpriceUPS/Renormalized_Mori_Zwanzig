function C = Ck_fill(irange,jrange,lrange,k,convo)
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
%  lrange  =  the range of k_z values (positive, negative, or zero)
%
%       k  =  the array of wavevectors (2Mx2Mx2Mx3)
%
%   convo  =  the convolution sum ( sum_{p+q=k} v_q*w_p^T )
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  C  =  the convolution C_k(v,w) for each Fourier mode in the range
%            [irange,jrange,lrange,:]

% find the size of the arrays we are using
s = size(convo);


t1 = zeros(s(1:end-1));
k1 = repmat(k(irange,jrange,lrange,1),[1,1,1,3]);
k2 = repmat(k(irange,jrange,lrange,2),[1,1,1,3]);
k3 = repmat(k(irange,jrange,lrange,3),[1,1,1,3]);
t1(irange,jrange,lrange,:) = -1j*(k1.*convo(irange,jrange,lrange,:,1) + k2.*convo(irange,jrange,lrange,:,2) + k3.*convo(irange,jrange,lrange,:,3));


t2 = zeros(s(1:end-1));

wave_vec_mag = sum(k(irange,jrange,lrange,:).^2,4);
wave_vec_mag(wave_vec_mag==0)=1;
wave_vec_mag = repmat(wave_vec_mag,[1,1,1,3]);

t1_1 = repmat(t1(irange,jrange,lrange,1),[1,1,1,3]);
t1_2 = repmat(t1(irange,jrange,lrange,2),[1,1,1,3]);
t1_3 = repmat(t1(irange,jrange,lrange,3),[1,1,1,3]);

B = k1.*t1_1 + k2.*t1_2 + k3.*t1_3;
t2(irange,jrange,lrange,:) = k(irange,jrange,lrange,:).*B./wave_vec_mag;




% retain only relevant terms
C = t1(irange,jrange,lrange,:) - t2(irange,jrange,lrange,:);
