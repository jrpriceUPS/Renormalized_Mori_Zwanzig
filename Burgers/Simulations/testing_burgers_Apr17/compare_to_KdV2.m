%clear all;close all;


addpath ../../nonlinear
addpath ../../simulation_functions
addpath ../../analysis

N = 64;
M = 3*N;
alpha = 1;
F_modes = [1:N,2*N:4*N+2,5*N+2:6*N];
G_modes = N+1:5*N+1;
k = [0:3*N-1,-3*N:-1].';
epsilon = 0;
degree = 4;

u = rand(N,1) + 1j*rand(N,1);

%compute Markov term
[t0,t0hat,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
t0 = t0(1:N);
[nonlin0,u_full2] = markov_term_KdV(u,M,N,alpha);


%compute t-model term
[t1,~,~] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes);
t1 = t1(1:N);
[nonlin1,uu_star] = tmodel_term_KdV(u_full2,nonlin0,alpha,F_modes);



%compute t^2-model term
[t2,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde] = t2model_term_Burgers(u_full,alpha,t0hat,t0tilde,F_modes,G_modes);
t2 = t2(1:N);
[nonlin2,uk3,uu,A,A_star,B,B_star,C,C_star,D,D_star] = t2model_term_complete_KdV(u_full2,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);




[t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde);
t3 = t3(1:N);
[nonlin3,uk6,E,E_star,F,F_star] = t3model_term_complete_KdV(alpha,F_modes,G_modes,k,epsilon,u_full2,uu,uu_star,uk3,A,A_star,B,B_star,C,C_star,D_star);




t4 = t4model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde,Ehat,Etilde,Fhat,Ftilde);
t4 = t4(1:N);
nonlin4 = t4model_term_complete_KdV(alpha,F_modes,G_modes,k,epsilon,u_full2,uu,uu_star,uk3,uk6,A,A_star,B,B_star,C,C_star,D,D_star,E,E_star,F,F_star);




nonlin0 = nonlin0(1:N);
err0 = sum(abs(t0 - nonlin0))./sum(abs(nonlin0))

nonlin1 = nonlin1(1:N);
err1 = sum(abs(t1 - nonlin1))./sum(abs(nonlin1))

nonlin2 = nonlin2(1:N);
err2 = sum(abs(t2 - nonlin2))./sum(abs(nonlin2))

nonlin3 = nonlin3(1:N);
err3 = sum(abs(t3 - nonlin3))./sum(abs(nonlin3))

nonlin4 = nonlin4(1:N);
err4 = sum(abs(t4 - nonlin4))./sum(abs(nonlin4))