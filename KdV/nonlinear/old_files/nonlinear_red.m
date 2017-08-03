function nonlin=nonlinear_red(u,alpha,M,N,F_modes,t)
%
%Computes the nonlinear part of the right hand side of the t-model of the
%KdV equation based upon a "full" model with M positive modes (M>N)

%compute Markov term
[nonlin0,u_full] = markov_term(u,M,N,alpha);

%compute t-model term
[nonlin1,~] = tmodel_term(u_full,nonlin0,alpha,F_modes);

%compute nonlinear part of right hand side of t-model for KdV
nonlin = nonlin0(1:N) + t*nonlin1(1:N);