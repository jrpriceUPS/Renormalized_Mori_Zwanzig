function err = log_periodic_cost(current,coeffs,N_list)

alpha = current(1);
beta = current(2);
gamma = current(3);

err = sum((beta*N_list.^alpha.*cos(gamma*log(N_list))-coeffs).^2);