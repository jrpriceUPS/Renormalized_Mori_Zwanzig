function [out] = memory_sign_odd(t_array,mem_vector)

    ind = find(mem_vector > 0 & abs(mem_vector) > eps);
    out.initial_time = t_array(ind(1));
    out.final_time = t_array(ind(end));
    out.value_inital = mem_vector(ind(1));
    out.value_final = mem_vector(ind(end));
    out.max_abs = max(abs(mem_vector(ind)));