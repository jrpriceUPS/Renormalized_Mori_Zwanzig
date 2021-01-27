function [c1_data,c2_data,c3_data,c4_data] = renormalize(alpha,N_list,u_list,t_list,tau_list)
%
%
% [c1_data,c2_data,c3_data,c4_data] = renormalize(alpha,N_list,u_list,t_list,tau_list)
%
% Computes optimal renormalization coefficients and tau for ROMs of degree
% n = 1,2,3,4 using the resolved data from "full" spectral solution
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%             alpha  =  degree of dispersion
%
%            N_list  =  set of resolutiosn for which to compute renormalization
%                       coefficients
%
%            u_list  =  Mxlength(t_list) array of the exact solution
%                       This is typically trimmed based on the resolved
%                       criterion prior to renormalization
%
%            t_list  =  the times associated with the u_list solutions
%                       This is typically trimmed based on the resolved
%                       criterion prior to renormalization
%
%           tau_list = degree of non-linearity of time dependence of
%                      renormalization coefficients
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
% 
% n = 1, 2, 3, 4
%
%         cn_data.cn = dictionary of n x length(N_list) set of coefficients for
%                      the t^n-model for each tau in tau_list
%
%         cn_data.En = dictionary of 1 x length(N_list) set of least squares 
%                      error for the coefficients for t^n-model for each tau in 
%                      tau_list
%   
%    cn_data.Axcondn = dictionary of 1 x length(N_list) set of condition
%                      number for A matrix for the coefficients for the t^n-model 
%                      for each tau in tau_list
%                  
%      cn_data.cn_op = n x length(N_list) set of optimal coefficients for the 
%                      t^n-model
%
%      cn_data.En_op = 1 x length(N_list) set of least squares error for the 
%                      optimal coefficients for t^n-model
%
% cn_data.Axcondn_op = 1 x length(N_list) set of condition number for A matrix 
%                      for the optimal coefficients for t^n-model
%
    
% Preallocate dictionaries
c1{1,length(N_list)} = [];
c2{1,length(N_list)} = [];
c3{1,length(N_list)} = [];
c4{1,length(N_list)} = [];

E1{1,length(N_list)} = [];
E2{1,length(N_list)} = [];
E3{1,length(N_list)} = [];
E4{1,length(N_list)} = [];

Axcond1{1,length(N_list)} = [];
Axcond2{1,length(N_list)} = [];
Axcond3{1,length(N_list)} = [];
Axcond4{1,length(N_list)} = [];

A4{1,length(N_list)} = [];

c1_op = zeros(2,length(N_list));
E1_op = zeros(1,length(N_list));
Axcond1_op = zeros(1,length(N_list));

c2_op = zeros(3,length(N_list));
E2_op = zeros(1,length(N_list));
Axcond2_op = zeros(1,length(N_list));

c3_op = zeros(4,length(N_list));
E3_op = zeros(1,length(N_list));
Axcond3_op = zeros(1,length(N_list));

c4_op = zeros(5,length(N_list));
E4_op = zeros(1,length(N_list));
Axcond4_op = zeros(1,length(N_list));
A4_op = zeros(4,4,length(N_list));

s = size(u_list);
N_full = s(1);
M_full = 3*N_full;
exact_derivative = zeros(size(u_list));

for i = 1:length(t_list)
    
    disp(sprintf('Calculating the exact derivative, the current time is t = %i',t_list(i)))
    u = u_list(1:N_full,i);
    [t0,~,~,~] = markov_term_Burgers(u,M_full,N_full,alpha,[],[]);
    exact_derivative(:,i) = t0(1:N_full).*conj(u) + conj(t0(1:N_full)).*u;
      
end

for i = 1:length(N_list)
    
    N = N_list(i);
    
    % Wave numbers
    F_modes = [1:N,2*N+1:4*N+1,5*N+1:6*N];
    G_modes = [N+1:5*N+1];
    M = 3*N;
    
    % Preallocated arrays in dictionary
    c1{i} = zeros(1,length(tau_list));
    c2{i} = zeros(2,length(tau_list));
    c3{i} = zeros(3,length(tau_list));
    c4{i} = zeros(4,length(tau_list));

    E1{i} = zeros(1,length(tau_list));
    E2{i} = zeros(1,length(tau_list));
    E3{i} = zeros(1,length(tau_list));
    E4{i} = zeros(1,length(tau_list));

    Axcond1{i} = zeros(1,length(tau_list));
    Axcond2{i} = zeros(1,length(tau_list));
    Axcond3{i} = zeros(1,length(tau_list));
    Axcond4{i} = zeros(1,length(tau_list));
    
    % Preallocated energy arrays
    markov_energy = zeros(N,length(t_list));
    tmodel_energy = zeros(N,length(t_list));
    t2model_energy = zeros(N,length(t_list));
    t3model_energy = zeros(N,length(t_list));
    t4model_energy = zeros(N,length(t_list));
    
    for j = 1:length(t_list)
        
        disp(sprintf('Calculating energy for N = %i, time t = %i',N,t_list(j)))
        u = u_list(1:N,j);
        
        %compute Markov term
        [t0,t0hat,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
        markov_energy(1:N,j) = t0(1:N).*conj(u) + conj(t0(1:N)).*u;
         
        %compute t-model term
        [t1,~,~] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes);
        tmodel_energy(1:N,j) = t1(1:N).*conj(u) + conj(t1(1:N)).*u;
        
        %compute t^2-model term
        [t2,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde] = t2model_term_Burgers(u_full,alpha,t0hat,t0tilde,F_modes,G_modes);
        t2model_energy(1:N,j) = t2(1:N).*conj(u) + conj(t2(1:N)).*u;
        
        %compute t^3-model term
        [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde);
        t3model_energy(1:N,j) = t3(1:N).*conj(u) + conj(t3(1:N)).*u;
       
        %compute t^4-model term
        t4 = t4model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde,Ehat,Etilde,Fhat,Ftilde);
        t4model_energy(1:N,j) = t4(1:N).*conj(u) + conj(t4(1:N)).*u;
      
    end
    
    for h = 1:length(tau_list)
        
        tau = tau_list(h);
        disp(sprintf('Currently calculating coefficients for N = %i: tau = %s',N,num2str(tau)))

        exact = exact_derivative(1:N,:);
        R0 = markov_energy(1:N,:);
        R1 = t_list.^(1-1*tau).*tmodel_energy(1:N,:);
        R2 = t_list.^(2-2*tau).*t2model_energy(1:N,:);
        R3 = t_list.^(3-3*tau).*t3model_energy(1:N,:);
        R4 = t_list.^(4-4*tau).*t4model_energy(1:N,:);

        % compute the RHS for the least squares solve
        RHS = R0 - exact;
        
        b = -sum(RHS(:).*R1(:));
        
        % construct the matrix for the least squares solve
        A11 = sum(R1(:).*R1(:));
        
        c1{i}(1,h) = A11\b;
        E1{i}(:,h) = sum(sum((R0+c1{i}(1,h)*R1-exact).^2,2),1);
        Axcond1{i}(1,h) = cond(A11);
        
        b = [b
            -sum(RHS(:).*R2(:))];
        
        % construct the matrix for the least squares solve
        A12 = sum(R1(:).*R2(:));
        
        A21 = A12;
        A22 = sum(R2(:).*R2(:));
        
        % solve the system and store the result
        A = [A11 A12
            A21 A22];
        
        c2{i}(1:2,h) = A\b;
        E2{i}(:,h) = sum(sum((R0+c2{i}(1,h)*R1+c2{i}(2,h)*R2-exact).^2,2),1);
        Axcond2{i}(1,h) = cond(A);
        
        b = [b
            -sum(RHS(:).*R3(:))];
        
        % construct the matrix for the least squares solve
        A13 = sum(R1(:).*R3(:));
        A23 = sum(R2(:).*R3(:));
        
        A31 = A13;
        A32 = A23;
        A33 = sum(R3(:).*R3(:));
        
        % solve the system and store the result
        A = [A11 A12 A13
            A21 A22 A23
            A31 A32 A33];
        
        c3{i}(1:3,h) = A\b;
        E3{i}(:,h) = sum(sum((R0+c3{i}(1,h)*R1+c3{i}(2,h)*R2+c3{i}(3,h)*R3-exact).^2,2),1);
        Axcond3{i}(1,h) = cond(A);
        
        b = [b
            -sum(RHS(:).*R4(:))];
        
        % construct the matrix for the least squares solve
        A14 = sum(R1(:).*R4(:));
        A24 = sum(R2(:).*R4(:));
        A34 = sum(R3(:).*R4(:));
        
        A41 = A14;
        A42 = A24;
        A43 = A34;
        A44 = sum(R4(:).*R4(:));
        
        % solve the system and store the result
        A = [A11 A12 A13 A14
            A21 A22 A23 A24
            A31 A32 A33 A34
            A41 A42 A43 A44];
        
        c4{i}(1:4,h) = A\b;
        E4{i}(:,h) = sum(sum((R0+c4{i}(1,h)*R1+c4{i}(2,h)*R2+c4{i}(3,h)*R3+c4{i}(4,h)*R4-exact).^2,2),1);
        Axcond4{i}(1,h) = cond(A);
        A4{i}(1:4,1:4,h) = A;
                
    end
   
    % Find the tau & renormalization coefficients corresponding to the minimum 
    % error for each N
    [~,ind1] = min(E1{1,i}(:));
    c1_op(1,i) = c1{1,i}(1,ind1);
    c1_op(2,i) = tau_list(ind1);
    E1_op(i) = E1{1,i}(ind1);
    Axcond1_op(i) = Axcond1{1,i}(ind1);
    
    [~,ind2] = min(E2{1,i}(:));
    c2_op(1:2,i) = c2{1,i}(1:2,ind2);
    c2_op(3,i) = tau_list(ind2);
    E2_op(i) = E2{1,i}(ind2);
    Axcond2_op(i) = Axcond2{1,i}(ind2);
    
    [~,ind3] = min(E3{1,i}(:));
    c3_op(1:3,i) = c3{1,i}(1:3,ind3);
    c3_op(4,i) = tau_list(ind3);
    E3_op(i) = E3{1,i}(ind3);
    Axcond3_op(i) = Axcond3{1,i}(ind3);
    
    [~,ind4] = min(E4{1,i}(:));
    c4_op(1:4,i) = c4{1,i}(1:4,ind4);
    c4_op(5,i) = tau_list(ind4);
    E4_op(i) = E4{1,i}(ind4);
    Axcond4_op(i) = Axcond4{1,i}(ind4);
    A4_op(1:4,1:4,i) = A4{1,i}(1:4,1:4,ind4);
    
end

% Build a structure to store all of the variables for each order
c1_data.c1 = c1;
c1_data.E1 = E1;
c1_data.Axcond1 = Axcond1;
c1_data.c1_op = c1_op;
c1_data.E1_op = E1_op;
c1_data.Axcond1_op = Axcond1_op;

c2_data.c2 = c2;
c2_data.E2 = E2;
c2_data.Axcond2 = Axcond2;
c2_data.c2_op = c2_op;
c2_data.E2_op = E2_op;
c2_data.Axcond2_op = Axcond2_op;

c3_data.c3 = c3;
c3_data.E3 = E3;
c3_data.Axcond3 = Axcond3;
c3_data.c3_op = c3_op;
c3_data.E3_op = E3_op;
c3_data.Axcond3_op = Axcond3_op;

c4_data.c4 = c4;
c4_data.E4 = E4;
c4_data.Axcond4 = Axcond4;
c4_data.c4_op = c4_op;
c4_data.E4_op = E4_op;
c4_data.Axcond4_op = Axcond4_op;
c4_data.A4_op = A4_op;




