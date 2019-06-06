function [dE,dE_CMA] = CMA_fixed_sign_check(u,N_list,t_list)
%
% A function to compute Delta E_1 through Delta E_4 to check whether the
% total flow in / out of the resolved modes has fixed sign for the CMA of
% 3D Euler
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       u  =  an MxMxMx3x4xlength(t_list) array of fully resolved Fourier modes
%
%  N_list  =  the resolutions to check
%
%  t_list  =  a vector of the times corresponding to the solution u
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%      dE  =  a length(t_list)x1 of the actual rate of energy entering /
%            leaving the resolved modes based on a "full" reference solution
%
%  dE_CMA  =  a length(t_list)x4 array of the rate of change of energy in the
%             resolved modes due to the t-model through t4-model

% initialize the output
dE = zeros(length(N_list),length(t_list),1);
dE_CMA = zeros(length(N_list),length(t_list),4);

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
    
    % construct output arrays of the proper size
    exact = zeros(2*N,2*N,2*N,3,t);
    R0 = zeros(2*N,2*N,2*N,3,t);
    R1 = zeros(2*N,2*N,2*N,3,t);
    R2 = zeros(2*N,2*N,2*N,3,t);
    R3 = zeros(2*N,2*N,2*N,3,t);
    R4 = zeros(2*N,2*N,2*N,3,t);
    
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
        R1(:,:,:,:,i) = t1_energy(res_exact,res_exact,res_exact,:);
        
        
        [t2,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2);
        t2 = u_squishify(t2,N);
        % compute the energy derivative due to the R_2 term
        t2_energy = t2.*conj(u_current) + conj(t2).*u_current;
        t2_energy = u_fullify(t2_energy,M_full);
        R2(:,:,:,:,i) = t2_energy(res_exact,res_exact,res_exact,:);
        
        
        
        [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde,a_tilde2);
        t3 = u_squishify(t3,N);
        % compute the energy derivative due to the R_3 term
        t3_energy = t3.*conj(u_current) + conj(t3).*u_current;
        t3_energy = u_fullify(t3_energy,M_full);
        R3(:,:,:,:,i) = t3_energy(res_exact,res_exact,res_exact,:);
        
        
        
        t4 = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde,a_tilde2);
        t4 = u_squishify(t4,N);
        % compute the energy derivative due to the R_4 term
        t4_energy = t4.*conj(u_current) + conj(t4).*u_current;
        t4_energy = u_fullify(t4_energy,M_full);
        R4(:,:,:,:,i) = t4_energy(res_exact,res_exact,res_exact,:);
        
        
        % record energy flow in reference solution
        exact_unspool = exact(:,:,:,:,i);
        dE(j,i) = sum(exact_unspool(:));
        
        % record total energy flow due to each term
        R1_unspool = R1(:,:,:,:,i);
        R2_unspool = R2(:,:,:,:,i);
        R3_unspool = R3(:,:,:,:,i);
        R4_unspool = R4(:,:,:,:,i);
        
        dE_CMA(j,i,1) = sum(R1_unspool(:));
        dE_CMA(j,i,2) = sum(R2_unspool(:));
        dE_CMA(j,i,3) = sum(R3_unspool(:));
        dE_CMA(j,i,4) = sum(R4_unspool(:));
        
        
        
    end
    
    figure(1)
    plot(t_list,dE(j,:))
    title(sprintf('Reference drain, N = %i',N),'fontsize',16)
    xlabel('time')
    ylabel('Resolved energy derivative')
    saveas(gcf,sprintf('total_drain%i',N),'png')
    close
    
    figure(2)
    for i = 1:4
        subplot(2,2,i)
        plot(t_list,dE_CMA(j,:,i))
        xlabel('time')
        ylabel('Resolved energy derivative')
        title(sprintf('t^%i-model drain, N = %i',i,N),'fontsize',16)
    end
    saveas(gcf,sprintf('CMA_drain%i',N),'png')
    close
    
end
