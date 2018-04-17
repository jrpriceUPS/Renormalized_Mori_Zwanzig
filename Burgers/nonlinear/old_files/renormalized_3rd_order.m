function nonlin=renormalized_3rd_order(u,t,simulation_params)
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
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
if coeffs(1) == 0
    nonlin1 = zeros(N,1);
else
    [nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);
end

%compute t^2-model term
if coeffs(1) == 0
    nonlin2 = zeros(N,1);
else
    [nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);
end

%compute t^3-model term
if coeffs(1) == 0
    nonlin3 = zeros(N,1);
else
    nonlin3 = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);
end

%compute nonlinear part of right hand side of t-model for KdV
nonlin = nonlin0(1:N) + coeffs(1)*t*nonlin1(1:N) - coeffs(2)*t^2/2*nonlin2(1:N) + coeffs(3)*t^3/6*nonlin3(1:N);