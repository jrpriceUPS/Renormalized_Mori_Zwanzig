function scaling_laws = renormalize_exponent(u,N_list,t_list,x0)
%
% A function to take a fully resolved u to compute optimal renormalization
% coefficients, including accounting for time exponent
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
%      x0  =  initial guess of coefficients for optimization
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  scaling_laws  =  a 4x3 array of the best fit scaling laws for these
%                   data for each ROM size and each term in the ROM

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

exact = cell(length(N_list),1);
R0 = cell(length(N_list),1);
R1 = cell(length(N_list),1);
R2 = cell(length(N_list),1);
R3 = cell(length(N_list),1);
R4 = cell(length(N_list),1);

N_array = cell(length(N_list),1);
t_array = cell(length(N_list),1);

for i = 1:length(t_list)
    for j = 1:length(N_list)
        N = N_list(j);
        t_array{j}(:,:,:,:,i) = repmat(t_list(i),2*N,2*N,2*N,3);
    end
end



% compute contribution to the derivative of the energy in each mode for
% each ROM at each resolution
for j = 1:length(N_list);
    
    % compute N and M for the current ROM
    N = N_list(j);
    M = 3*N;
    exact{j} = zeros(2*N,2*N,2*N,3,t);
    R0{j} = zeros(2*N,2*N,2*N,3,t);
    R1{j} = zeros(2*N,2*N,2*N,3,t);
    R2{j} = zeros(2*N,2*N,2*N,3,t);
    R3{j} = zeros(2*N,2*N,2*N,3,t);
    R4{j} = zeros(2*N,2*N,2*N,3,t);
    
    N_array{j} = repmat(N,2*N,2*N,2*N,3,t);
    
    
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
        disp(sprintf('Current time is t = %i out of %i\n',current_t,t_list(end)))
        % find the exact derivative for the terms in this ROM
        exact_full = u_fullify(squeeze(exact_everything(:,:,:,:,:,i)),M_full);
        exact{j}(:,:,:,:,i) = exact_full(res_exact,res_exact,res_exact,:);
        
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
        R0{j}(:,:,:,:,i) = t0_energy(res_exact,res_exact,res_exact,:);
        
        [~,t1hat,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde,a_tilde2);
        t1 = u_squishify(t1hat,N);
        
        % compute the energy derivative due to the R_1 term
        t1_energy = t1.*conj(u_current) + conj(t1).*u_current;
        t1_energy = u_fullify(t1_energy,M_full);
        R1{j}(:,:,:,:,i) = t1_energy(res_exact,res_exact,res_exact,:);
        
        
        [t2,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2);
        t2 = u_squishify(t2,N);
        
        % compute the energy derivative due to the R_2 term
        t2_energy = t2.*conj(u_current) + conj(t2).*u_current;
        t2_energy = u_fullify(t2_energy,M_full);
        R2{j}(:,:,:,:,i) = t2_energy(res_exact,res_exact,res_exact,:);
        
        
        
        [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde,a_tilde2);
        t3 = u_squishify(t3,N);
        
        % compute the energy derivative due to the R_3 term
        t3_energy = t3.*conj(u_current) + conj(t3).*u_current;
        t3_energy = u_fullify(t3_energy,M_full);
        R3{j}(:,:,:,:,i) = t3_energy(res_exact,res_exact,res_exact,:);
        
        
        
        t4 = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde,a_tilde2);
        t4 = u_squishify(t4,N);
        
        % compute the energy derivative due to the R_4 term
        t4_energy = t4.*conj(u_current) + conj(t4).*u_current;
        t4_energy = u_fullify(t4_energy,M_full);
        R4{j}(:,:,:,:,i) = t4_energy(res_exact,res_exact,res_exact,:);
    end
end
    
scaling_laws = fminsearch(@(x)  optimize(x,N_array,t_array,exact,R0,R1,R2,R3,R4),x0);