function [scaling_laws,corr_coeff] = create_scaling_laws(N_list,c)
%
%  [scaling_laws,corr_coeff] = create_scaling_laws(N_list,c)
%
%  Computes scaling laws and correlation coefficients from optimal
%  renormalization coefficients
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%  N_list  =  list of resolutions coefficients correspond to
%
%       c  =  set of optimal coefficients
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  scaling_laws  =  set of scaling laws associated with each degree of ROM
%
%    corr_coeff  =  goodness of fit for each variable

% check the size of the coefficient array
s = size(c);

% if the length is 3, we have coefficients calculated as a function of time
if length(s)==3
    corr_coeff = zeros(s(1),s(3));
    
    scaling_laws = zeros(s(1),2,s(3));
    
    % compute the scaling law fit for all degrees of coefficients for all
    % timesteps
    for i = 1:s(1)
        for j = 1:s(3)
            
            r = polyfit(log(N_list),log(abs(squeeze(c(i,:,j)))),1);
            scaling_laws(i,:,j) = r;
            corr_mat = corrcoef(log(N_list),log(squeeze(c(i,:,j))));
            corr_coeff(i,j) = corr_mat(1,2)^2;
            
        end
    end
end

% if the length is 2, we have only one set of coefficients for each
% resolution and degree
if length(s) == 2
    corr_coeff = zeros(s(1),1);
    
    % compute the scaling law fit for all degrees of coefficients
    for i = 1:s(1)
        
        r = polyfit(log(N_list),log(abs(squeeze(c(i,:)))),1);
        scaling_laws(i,:) = r;
        corr_mat = corrcoef(log(N_list),log(squeeze(c(i,:))));
        corr_coeff(i) = corr_mat(1,2)^2;
        
        
    end
end