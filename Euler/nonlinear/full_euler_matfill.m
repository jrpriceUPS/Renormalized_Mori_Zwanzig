function [du_dt,term1,term2] = full_euler_matfill(irange,jrange,lrange,k,convo)
%
% A function to compute the derivative of a specified range of wavenumbers
% for the full model
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
%   convo  =  the convolution sum ( sum_{p+q=k} u_q*u_p^T )
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  du_dt  =  the time derivative of each Fourier mode in the range
%            [irange,jrange,lrange,:]
%
%  term1  =  the term -i * sum_{p+q = k} u_q*u_p^T * k
%
%  term2  =  the term i * k / |k|^2 * sum_{p+q = k} (k dot u_p)(k dot u_q)

% find the size of the arrays we are using
s = size(convo);

% construct full output arrays (we'll only update specific entries)
term1 = zeros(s(1:end-1));
term2 = zeros(s(1:end-1));
du_dt = zeros(s(1:end-1));

% loop through the entries we want to update
for i = irange
    for j = jrange
        for l = lrange
            
            if i == 1 && j == 1 && l == 1
                
                term1(i,j,l,:) = zeros(3,1);
                term2(i,j,l,:) = zeros(3,1);
                du_dt(i,j,l,:) = term1(i,j,l,:) - term2(i,j,l,:);
                
            else
                
                wave_vec = squeeze(k(i,j,l,:)); % current wavevector k
                current_mat = squeeze(convo(i,j,l,:,:)); % associated matrix
                
                term1(i,j,l,:) = -1j*current_mat*wave_vec;
                term2(i,j,l,:) = wave_vec*wave_vec.'*squeeze(term1(i,j,l,:))/norm(wave_vec)^2;
                du_dt(i,j,l,:) = term1(i,j,l,:) - term2(i,j,l,:);
                
            end
        end
    end
end

% retain only relevant terms
term1 = term1(irange,jrange,lrange,:);
term2 = term2(irange,jrange,lrange,:);
du_dt = du_dt(irange,jrange,lrange,:);