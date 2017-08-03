function nonlin=multi_renormalized_3rd_order(u,t,simulation_params)
%
%Computes the nonlinear part of the right hand side of the t^3-model of the
%KdV equation based upon a "full" model with M positive modes (M>N)

%gather parameters needed for simulation
alpha = simulation_params.alpha;
F_modes = simulation_params.F_modes;
G_modes = simulation_params.G_modes;
k = simulation_params.k;
coeffs = simulation_params.coeffs;
N = simulation_params.N;
M = simulation_params.M;
epsilon =simulation_params.epsilon;


%compute Markov term
[nonlin0_nondispers,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
if coeffs(1) == 0
    nonlin1_nondispers = zeros(N,1);
else
    [nonlin1_nondispers,uu_star] = tmodel_term(u_full,nonlin0_nondispers,alpha,F_modes);
end

%compute t^2-model term
if coeffs(2) == 0
    nonlin2_nondispers = zeros(N,1);
else
    [nonlin2_nondispers,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0_nondispers,uu_star,alpha,F_modes,G_modes,k,epsilon);
end

%compute t^3-model term
if coeffs(3) == 0
    nonlin3_nondispers = zeros(N,1);
else
    nonlin3_nondispers = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);
end



%compute t^2-model term
if coeffs(2) == 0
    nonlin2_dispers = zeros(N,1);
else
    [nonlin2_dispers,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0_nondispers,uu_star,alpha,F_modes,G_modes,k,epsilon);
end

%compute t^3-model term
if coeffs(3) == 0
    nonlin3_dispers = zeros(N,1);
else
    nonlin3_dispers = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);
end

nonlin2_dispers = nonlin2_dispers - nonlin2_nondispers;
nonlin3_dispers = nonlin3_dispers - nonlin3_nondispers;

%compute nonlinear part of right hand side of t-model for KdV
nonlin = nonlin0_nondispers(1:N) + coeffs(1,1)*t*nonlin1_nondispers(1:N) - coeffs(2,1)*t^2/2*nonlin2_nondispers(1:N) + coeffs(3,1)*t^3/6*nonlin3_nondispers(1:N)...
    - coeffs(2,1)*t^2/2*nonlin2_dispers(1:N) + coeffs(3,2)*t^3/6*nonlin3_dispers(1:N);