function [scaling_laws,corr_coeff] = create_scaling_laws(N_list,c)

s = size(c);

if length(s)==3
    corr_coeff = zeros(s(1),s(3));
    
    scaling_laws = zeros(s(1),2,s(3));
    
    for i = 1:s(1)
        for j = 1:s(3)
            
            r = polyfit(log(N_list),log(abs(squeeze(c(i,:,j)))),1);
            scaling_laws(i,:,j) = r;
            corr_mat = corrcoef(log(N_list),log(squeeze(c(i,:,j))));
            corr_coeff(i,j) = corr_mat(1,2)^2;
            
        end
    end
end

if length(s) == 2
    corr_coeff = zeros(s(1),1);
    
    for i = 1:s(1)
        
        r = polyfit(log(N_list),log(abs(squeeze(c(i,:)))),1);
        scaling_laws(i,:) = r;
        corr_mat = corrcoef(log(N_list),log(squeeze(c(i,:))));
        corr_coeff(i) = corr_mat(1,2)^2;
        
        
    end
end