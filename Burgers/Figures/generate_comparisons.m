function [times_array,energies_array,error_array] = generate_comparisons(N_list,endtime)

times_array = cell(length(N_list),1);
energies_array = cell(length(N_list),1);
error_array = cell(length(N_list),1);

for i = 1:length(N_list)
    [times,energies,errors] = error_test(N_list(i),1,endtime);
    times_array{i} = times;
    energies_array{i} = energies;
    error_array{i} = errors;
end