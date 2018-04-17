function renormalization_params = multi_renormalized_3rd_order_energy_dispers(simulation_params,renormalization_params)

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
coeffs = simulation_params.coeffs(:,1);


%compute Markov term
[nonlin0_nondispers,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1_nondispers,uu_star] = tmodel_term(u_full,nonlin0_nondispers,alpha,F_modes);

%compute t^2-model term
[nonlin2_nondispers,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0_nondispers,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute t^3-model term
nonlin3_nondispers = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);


%compute nonlinear part of right hand side of t-model for KdV
nonlinear_markov = sum(nonlin0_nondispers(1:N).*conj(u) + conj(nonlin0_nondispers(1:N)).*u);
nonlinear_markov = nonlinear_markov + coeffs(1)*(t*sum(nonlin1_nondispers(1:N).*conj(u) + conj(nonlin1_nondispers(1:N)).*u));
nonlinear_markov = nonlinear_markov + coeffs(2)*(-t^2/2*sum(nonlin2_nondispers(1:N).*conj(u) + conj(nonlin2_nondispers(1:N)).*u));
nonlinear_markov = nonlinear_markov + coeffs(3)*(t^3/6*sum(nonlin3_nondispers(1:N).*conj(u) + conj(nonlin3_nondispers(1:N)).*u));

renormalization_params.markov_term(index) = nonlinear_markov;



epsilon = simulation_params.epsilon;


%compute Markov term
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[~,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);

%compute t^2-model term
[nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute t^3-model term
nonlin3 = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);


nonlin2 = nonlin2 - nonlin2_nondispers;
nonlin3 = nonlin3 - nonlin3_nondispers;

nonlinear_terms = zeros(1,3);
%compute nonlinear part of right hand side of t-model for KdV
nonlinear_terms(1) = -t^2/2*sum(nonlin2(1:N).*conj(u) + conj(nonlin2(1:N)).*u);
nonlinear_terms(2) = t^3/6*sum(nonlin3(1:N).*conj(u) + conj(nonlin3(1:N)).*u);

renormalization_params.nonlinear_terms(index,:) = nonlinear_terms;