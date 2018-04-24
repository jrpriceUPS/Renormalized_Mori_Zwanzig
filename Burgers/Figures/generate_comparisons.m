clear all;close all;

N_list = 4:2:24;

t_array = zeros(4,2,length(N_list));
for i = N_list
    t = coeff_test(i,1,10);
    t_array(:,:,i) = t;
end