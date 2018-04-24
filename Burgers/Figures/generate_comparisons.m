function error_array = generate_comparisons(N_list,endtime)

error_array = cell(length(N_list),4);

for i = 1:length(N_list)
    [errc1B,errc2B,errc3B,errc4B] = error_test(N_list(i),1,endtime);
    error_array{i,1} = errc1B;
    error_array{i,2} = errc2B;
    error_array{i,3} = errc3B;
    error_array{i,4} = errc4B;
end