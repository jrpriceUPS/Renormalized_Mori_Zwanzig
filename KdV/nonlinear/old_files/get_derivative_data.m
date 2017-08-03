function [nonlin0,nonlin1,nonlin2,nonlin3]=get_derivative_data(u,M,N,epsilon,alpha,F_modes,G_modes,k,t)
%
%Computes the nonlinear part of the right hand side of the t^3-model of the
%KdV equation based upon a "full" model with M positive modes (M>N)

%compute Markov term
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1,uu_star] = tmodel_term(M,N,u_full,M/N*nonlin0,alpha,F_modes);

%compute t^2-model term
[nonlin2,uk,uu,uuk_raw,uuk_star,uuu_raw,uuu_star] = t2model_term(M,N,u_full,M/N*nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);

%compute t^3-model term
nonlin3 = t3model_term(M,N,alpha,F_modes,G_modes,k,epsilon,u_full,uu,uu_star,uk,uuk_raw,uuk_star,uuu_raw,uuu_star);

%include the time dependence
nonlin1 = nonlin1*t;
nonlin2 = nonlin2*t^2;
nonlin3 = nonlin3*t^3;