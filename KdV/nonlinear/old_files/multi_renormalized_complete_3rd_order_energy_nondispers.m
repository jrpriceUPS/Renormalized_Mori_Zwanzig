function renormalization_params = multi_renormalized_complete_3rd_order_energy_nondispers(simulation_params,renormalization_params)

index = renormalization_params.index;
N = renormalization_params.N_list(index);
M = 4*N+1;
u = simulation_params.u(1:N);
epsilon = 0;
alpha = simulation_params.alpha;
F_modes = simulation_params.F_modes{index};
G_modes = simulation_params.G_modes{index};
k = simulation_params.k{index};
t = simulation_params.current_time;


%compute Markov term
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);

%compute t^2-model term
[nonlin2,uk,uu,A,A_star,B,B_star,C,C_star,D_star] = t2model_term_complete(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute t^3-model term
nonlin3 = t3model_term_complete(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,A,A_star,B,B_star,C,C_star,D_star);


nonlinear_terms = zeros(1,3);
%compute nonlinear part of right hand side of t-model for KdV
nonlinear_markov = sum(nonlin0(1:N).*conj(u) + conj(nonlin0(1:N)).*u);
nonlinear_terms(1) = t*sum(nonlin1(1:N).*conj(u) + conj(nonlin1(1:N)).*u);
nonlinear_terms(2) = -t^2/2*sum(nonlin2(1:N).*conj(u) + conj(nonlin2(1:N)).*u);
nonlinear_terms(3) = t^3/6*sum(nonlin3(1:N).*conj(u) + conj(nonlin3(1:N)).*u);

renormalization_params.markov_term(index) = nonlinear_markov;
renormalization_params.nonlinear_terms(index,:) = nonlinear_terms;