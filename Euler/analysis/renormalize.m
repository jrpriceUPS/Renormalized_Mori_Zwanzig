function [c1_data,c2_data,c3_data,c4_data] = renormalize(u,N_list,t_list,tau_list)
%
% A function to take a fully resolved u to compute optimal renormalization
% coefficients for several ROMs
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       u  =  an MxMxMx3x4xlength(t_list) array of fully resolved Fourier modes
%
%  N_list  =  a vector of resolutions (less than M) for which to compute
%             optimal renormalization coefficients
%
%  t_list  =  a vector of the times corresponding to the solution u
%
%  tau_list = degree of non-linearity of time dependence of
%              renormalization coefficients
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%   n = 1, 2, 3, 4
%   cn_data.cn = dictionary of n x length(N_list) set of coefficients for
%   the t^n-model for each tau in tau_list
%   cn_data.En = dictionary of 1 x length(N_list) set of least squares error for the coefficients for t^n-model for each tau in tau_list
%   cn_data.Axcondn = dictionary of 1 x length(N_list) set of condition
%   number for A matrix for the coefficients for the t^n-model for each tau
%   in tau_list
%   cn_data.cn_op = n x length(N_list) set of optimal coefficients for the t^n-model
%   cn_data.En_op = 1 x length(N_list) set of least squares error for the optimal coefficients for t^n-model
%   cn_data.Axcondn_op = 1 x length(N_list) set of condition number for A matrix for the optimal coefficients for t^n-model
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

% extract data about the size of the full array
s = size(u);
N_full = s(1);
M_full = 3*N_full;
t = s(6);

% initialize an array into which the exact derivative of the energy in each
% mode will be stored
exact_everything = zeros(s);

% make k array
k_vec_full = [0:M_full-1,-M_full:1:-1];
[kx_full,ky_full,kz_full] = ndgrid(k_vec_full,k_vec_full,k_vec_full);
k_full = zeros(2*M_full,2*M_full,2*M_full,3);
k_full(:,:,:,1) = kx_full;
k_full(:,:,:,2) = ky_full;
k_full(:,:,:,3) = kz_full;

% make vectors of resolved and unresolved modes in the full simulation
a_full = 2:M_full;
b_full = 2*M_full:-1:M_full+2;

% compute exact derivative of the energy in each mode at each timestep
for i = 1:t
    
    % compute and display the current time
    disp('Currently analyzing exact solution')
    disp(sprintf('Current time is t = %i out of %i\n',t_list(i),t_list(end)))
    
    % extract the current u
    u_current = squeeze(u(:,:,:,:,:,i));
    
    % construct the full version of it
    u_full = u_fullify(u_current,M_full);
    
    % compute the time derivative and convert it into the same shape as u
    % itself
    du_dt = Ck(u_full,u_full,a_full,b_full,k_full,[],[]);
    du_dt = u_squishify(du_dt,N_full);
    
    % compute the energy derivative of each mode for this timestep
    exact_everything(:,:,:,:,:,i) = du_dt.*conj(u_current) + conj(du_dt).*u_current;
end

% compute contribution to the derivative of the energy in each mode for
% each ROM at each resolution
for j = 1:length(N_list)
    
    % compute N and M for the current ROM
    N = N_list(j);
    M = 3*N;
    
    % Preallocated arrays in dictionary
    c1{j} = zeros(1,length(tau_list));
    c2{j} = zeros(2,length(tau_list));
    c3{j} = zeros(3,length(tau_list));
    c4{j} = zeros(4,length(tau_list));

    E1{j} = zeros(1,length(tau_list));
    E2{j} = zeros(1,length(tau_list));
    E3{j} = zeros(1,length(tau_list));
    E4{j} = zeros(1,length(tau_list));

    Axcond1{j} = zeros(1,length(tau_list));
    Axcond2{j} = zeros(1,length(tau_list));
    Axcond3{j} = zeros(1,length(tau_list));
    Axcond4{j} = zeros(1,length(tau_list));
    
    % construct output arrays of the proper size
    exact = zeros(2*N,2*N,2*N,3,t);
    R0 = zeros(2*N,2*N,2*N,3,t);
    R1_nt = zeros(2*N,2*N,2*N,3,t);
    R2_nt = zeros(2*N,2*N,2*N,3,t);
    R3_nt = zeros(2*N,2*N,2*N,3,t);
    R4_nt = zeros(2*N,2*N,2*N,3,t);
    t_array = zeros(2*N,2*N,2*N,3);
    
    % compute the indices corresponding to terms in this ROM from the exact
    % data computed above
    res_exact = [1:N,2*M_full - N+1:2*M_full];
    
    % make k array for the ROM
    k_vec = [0:M-1,-M:1:-1];
    [kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
    k = zeros(2*M,2*M,2*M,3);
    k(:,:,:,1) = kx;
    k(:,:,:,2) = ky;
    k(:,:,:,3) = kz;
    
    % compute indices of resolved and unresolved modes
    a = 2:M;
    b = 2*M:-1:M+2;
    a_tilde = N+1:M;
    a_tilde2 = 2*N+1:M;
    
    % loop through every timestep
    for i = 1:t
        disp(sprintf('Current resolution is N = %i',N))

        current_t = t_list(i);
        t_array(:,:,:,:,i) = repmat(t_list(i),2*N,2*N,2*N,3);

        disp(sprintf('Current time is t = %i out of %i\n',current_t,t_list(end)))
        
        % find the exact derivative for the terms in this ROM
        exact_full = u_fullify(squeeze(exact_everything(:,:,:,:,:,i)),M_full);
        exact(:,:,:,:,i) = exact_full(res_exact,res_exact,res_exact,:);
       
        % extract u and convert it into a full array with the proper
        % elements eliminated
        temp_u = squeeze(u(:,:,:,:,:,i));
        temp_u_full = u_fullify(temp_u,M_full);
        u_current = u_squishify(temp_u_full,N);
        u_full = u_fullify(u_current,M);
        
        % compute the R_k^i terms for i = 1,2,3,4
        [~,t0hat,t0tilde] = markov_term(u_full,a,b,k,a_tilde,a_tilde2);
        t0 = u_squishify(t0hat,N);
        % compute the energy derivative due to the R_0 term
        t0_energy = t0.*conj(u_current) + conj(t0).*u_current;
        t0_energy = u_fullify(t0_energy,M_full);
        R0(:,:,:,:,i) = t0_energy(res_exact,res_exact,res_exact,:);
        
        [~,t1hat,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde,a_tilde2);
        t1 = u_squishify(t1hat,N);
        
        % compute the energy derivative due to the R_1 term
        t1_energy = t1.*conj(u_current) + conj(t1).*u_current;
        t1_energy = u_fullify(t1_energy,M_full);
        R1_nt(:,:,:,:,i) = t1_energy(res_exact,res_exact,res_exact,:);
        
        [t2hat,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2);
        t2 = u_squishify(t2hat,N);
        
        % compute the energy derivative due to the R_2 term
        t2_energy = t2.*conj(u_current) + conj(t2).*u_current;
        t2_energy = u_fullify(t2_energy,M_full);
        R2_nt(:,:,:,:,i) = t2_energy(res_exact,res_exact,res_exact,:);
        
        [t3hat,Ehat,Etilde,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde,a_tilde2);
        t3 = u_squishify(t3hat,N);
        
        % compute the energy derivative due to the R_3 term
        t3_energy = t3.*conj(u_current) + conj(t3).*u_current;
        t3_energy = u_fullify(t3_energy,M_full);
        R3_nt(:,:,:,:,i) = t3_energy(res_exact,res_exact,res_exact,:);
        
        t4hat = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde,a_tilde2);
        t4 = u_squishify(t4hat,N);
        
        % compute the energy derivative due to the R_4 term
        t4_energy = t4.*conj(u_current) + conj(t4).*u_current;
        t4_energy = u_fullify(t4_energy,M_full);
        R4_nt(:,:,:,:,i) = t4_energy(res_exact,res_exact,res_exact,:);
        
    end
    
    for h = 1:length(tau_list)
        
        tau = tau_list(h);
        disp(sprintf('Currently calculating coefficients for N = %i: tau = %s',N,num2str(tau)))
        
        R1 = t_array.^(1-1*tau).*R1_nt;
        R2 = t_array.^(2-2*tau).*R2_nt;
        R3 = t_array.^(3-3*tau).*R3_nt;
        R4 = t_array.^(4-4*tau).*R4_nt;

        % compute the RHS for the least squares solve
        RHS = R0 - exact;
        
        b = -sum(RHS(:).*R1(:));
        
        % construct the matrix for the least squares solve
        A11 = sum(R1(:).*R1(:));
        
        c1{j}(1,h) = A11\b;
        E1{j}(:,h) = sum((R0(:)+c1{j}(1,h).*R1(:)-exact(:)).^2);
        Axcond1{j}(1,h) = cond(A11);
        
        b = [b
            -sum(RHS(:).*R2(:))];
        
        % construct the matrix for the least squares solve
        A12 = sum(R1(:).*R2(:));
        
        A21 = A12;
        A22 = sum(R2(:).*R2(:));
        
        % solve the system and store the result
        A = [A11 A12
            A21 A22];
        
        c2{j}(1:2,h) = A\b;
        E2{j}(:,h) = sum((R0(:)+c2{j}(1,h).*R1(:)+c2{j}(2,h).*R2(:)-exact(:)).^2);
        Axcond2{j}(1,h) = cond(A);
               
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
        
        c3{j}(1:3,h) = A\b;
        E3{j}(:,h) = sum((R0(:)+c3{j}(1,h).*R1(:)+c3{j}(2,h).*R2(:)+c3{j}(3,h).*R3(:)-exact(:)).^2);
        Axcond3{j}(1,h) = cond(A);
        
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
        
        c4{j}(1:4,h) = A\b;
        E4{j}(:,h) = sum((R0(:)+c4{j}(1,h).*R1(:)+c4{j}(2,h).*R2(:)+c4{j}(3,h).*R3(:)+c4{j}(4,h).*R4(:)-exact(:)).^2);
        Axcond4{j}(1,h) = cond(A);
        
    end
    
    % Find the tau & renormalization coefficients corresponding to the minimum
    % error for each N
    [~,ind1] = min(E1{1,j}(:));
    c1_op(1,j) = c1{1,j}(1,ind1);
    c1_op(2,j) = tau_list(ind1);
    E1_op(j) = E1{1,j}(ind1);
    Axcond1_op(j) = Axcond1{1,j}(ind1);
    
    [~,ind2] = min(E2{1,j}(:));
    c2_op(1:2,j) = c2{1,j}(1:2,ind2);
    c2_op(3,j) = tau_list(ind2);
    E2_op(j) = E2{1,j}(ind2);
    Axcond2_op(j) = Axcond2{1,j}(ind2);
    
    [~,ind3] = min(E3{1,j}(:));
    c3_op(1:3,j) = c3{1,j}(1:3,ind3);
    c3_op(4,j) = tau_list(ind3);
    E3_op(j) = E3{1,j}(ind3);
    Axcond3_op(j) = Axcond3{1,j}(ind3);
    
    [~,ind4] = min(E4{1,j}(:));
    c4_op(1:4,j) = c4{1,j}(1:4,ind4);
    c4_op(5,j) = tau_list(ind4);
    E4_op(j) = E4{1,j}(ind4);
    Axcond4_op(j) = Axcond4{1,j}(ind4);
    
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
    
