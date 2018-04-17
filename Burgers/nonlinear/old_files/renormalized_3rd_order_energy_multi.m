function renormalization_params=renormalized_3rd_order_energy_multi(simulation_params,renormalization_params)
%
%Computes the nonlinear part of the right hand side of the t^3-model of the
%KdV equation based upon a "full" model with M positive modes (M>N)

N_list = renormalization_params.N_list;
N = N_list(1);
M = 2*N+1;
u = simulation_params.u(1:N);
epsilon = simulation_params.epsilon;
alpha = simulation_params.alpha;
F_modes = simulation_params.F_modes{1};
G_modes = simulation_params.G_modes{1};
k = simulation_params.k{1};
t = simulation_params.current_time;


%compute Markov term
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);

%compute t^2-model term
[nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute t^3-model term
nonlin3 = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);

nonlinear_markov = zeros(length(N_list),1);
nonlinear_terms = zeros(length(N_list),3);

for i = 1:length(N_list)
    N0 = N_list(i);
    nonlinear_markov(i) = sum(nonlin0(1:N0).*conj(u(1:N0)) + conj(nonlin0(1:N0)).*u(1:N0));
    nonlinear_terms(i,1) = t*sum(nonlin1(1:N0).*conj(u(1:N0)) + conj(nonlin1(1:N0)).*u(1:N0));
    nonlinear_terms(i,2) = -t^2/2*sum(nonlin2(1:N0).*conj(u(1:N0)) + conj(nonlin2(1:N0)).*u(1:N0));
    nonlinear_terms(i,3) = t^3/6*sum(nonlin3(1:N0).*conj(u(1:N0)) + conj(nonlin3(1:N0)).*u(1:N0));
end

renormalization_params.markov_term = nonlinear_markov;
renormalization_params.nonlinear_terms = nonlinear_terms;