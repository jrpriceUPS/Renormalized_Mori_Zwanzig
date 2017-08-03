function [nonlinear_markov,nonlinear_terms]=renormalized_complete_3rd_order_energy_old(u,M,N,epsilon,alpha,F_modes,G_modes,k,t)
%
%Computes the nonlinear part of the right hand side of the t^3-model of the
%KdV equation based upon a "full" model with M positive modes (M>N)


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