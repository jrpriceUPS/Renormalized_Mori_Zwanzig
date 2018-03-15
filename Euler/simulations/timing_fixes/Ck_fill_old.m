function C = Ck_fill_old(irange,jrange,lrange,k,convo)
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

% construct full output arrays (we'll only update specific entries)
C = zeros(s(1:end-1));


% loop through the entries we want to update
for i = irange
    for j = jrange
        for l = lrange
            
            if i == 1 && j == 1 && l == 1
                
                term1 = zeros(3,1);
                term2 = zeros(3,1);
                C(i,j,l,:) = term1 - term2;
                
            else
                
                wave_vec = squeeze(k(i,j,l,:)); % current wavevector k
                current_mat = squeeze(convo(i,j,l,:,:)); % associated matrix
                
                term1 = -1j*current_mat*wave_vec;
                term2 = wave_vec*wave_vec.'*squeeze(term1)/norm(wave_vec)^2;
                C(i,j,l,:) = term1 - term2;
                
            end
        end
    end
end


% retain only relevant terms
C = C(irange,jrange,lrange,:);