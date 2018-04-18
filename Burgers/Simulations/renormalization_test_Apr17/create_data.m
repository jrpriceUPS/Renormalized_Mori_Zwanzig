function create_data(alpha,num_points,endtime,dt,howoften)

[t_list,u_list] = upwind_burgers(alpha,num_points,endtime,dt,howoften);
save u_list u_list
save t_list t_list

N = num_points/2;
F_modes = [1:N,2*N:4*N+2,5*N+2:6*N];
G_modes = N+1:5*N+1;
M = 3*N;
exact_derivative = zeros(size(u_list));

for i = 1:length(t_list)
    disp(sprintf('The current time is t = %i',t_list(i)))
    u = u_list(:,i);
    t0 = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
    exact_derivative(:,i) = t0(1:N).*conj(u) + conj(t0(1:N)).*u;
    
end

save exact_derivative exact_derivative