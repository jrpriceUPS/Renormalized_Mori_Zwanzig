function nonlin=nonlinear_3rd_order(u,M,N,epsilon,alpha,F_modes,G_modes,k,t)
%
%Computes the nonlinear part of the right hand side of the t^3-model of the
%KdV equation based upon a "full" model with M positive modes (M>N)


%compute Markov term
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);

%compute t^2-model term
[nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute t^3-model term
nonlin3 = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);

%compute nonlinear part of right hand side of t-model for KdV
nonlin = nonlin0(1:N) + t*nonlin1(1:N) - t^2/2*nonlin2(1:N) + t^3/6*nonlin3(1:N);