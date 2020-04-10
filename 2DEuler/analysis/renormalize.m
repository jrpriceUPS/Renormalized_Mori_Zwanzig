function [coeff_array,scaling_laws,r] = renormalize(u,N_list,t_list,time)
%
% A function to take a fully resolved u to compute optimal renormalization
% coefficients for several ROMs
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       u  =  an MxMx2x2xlength(t_list) array of fully resolved Fourier modes
%
%  N_list  =  a vector of resolutions (less than M) for which to compute
%             optimal renormalization coefficients
%
%  t_list  =  a vector of the times corresponding to the solution u
%
%    time  =  a logical variable (1 if algebraic time dependence of
%             renormalization coefficients is retained, 0 if they are not)
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
%  coeff_array   =  a 4xlength(N_list)x4 array of the optimal coefficients of
%                   each term in the ROM model for each resolution and ROM
%                   degree
%
%  scaling_laws  =  a 4x2x4 array of the best fit scaling laws for these
%                   data for each ROM size and each term in the ROM
%
%             r  =  a 4x4 array of the correlation coefficient of the best
%                   fit of each coefficient of each ROM term for each ROM size

% initialize the output
coeff_array = zeros(4,length(N_list),4);
r = zeros(4,4);

% extract data about the size of the full array
s = size(u);
N_full = s(1);
M_full = 3*N_full;
t = s(5);

% initialize an array into which the exact derivative of the energy in each
% mode will be stored
exact_everything = zeros(s);

% make k array
k_vec_full = [0:M_full-1,-M_full:1:-1];
[kx_full,ky_full] = ndgrid(k_vec_full,k_vec_full);
k_full = zeros(2*M_full,2*M_full,2);
k_full(:,:,1) = kx_full;
k_full(:,:,2) = ky_full;

% make vectors of resolved and unresolved modes in the full simulation
a_full = 2:M_full;
b_full = 2*M_full:-1:M_full+2;

Legend = cell(length(N_list),1);

% compute exact derivative of the energy in each mode at each timestep
for i = 1:t
    
    % compute and display the current time
    disp('Currently analyzing exact solution')
    fprintf('Current time is t = %i out of %i\n\n',t_list(i),t_list(end))
    
    % extract the current u
    u_current = squeeze(u(:,:,:,:,i));
    
    % construct the full version of it
    u_full = u_fullify(u_current,M_full);
    
    % compute the time derivative and convert it into the same shape as u
    % itself
    du_dt = Ck(u_full,u_full,a_full,b_full,k_full,[],[]);
    du_dt = u_squishify(du_dt,N_full);
    
    % compute the energy derivative of each mode for this timestep
    exact_everything(:,:,:,:,i) = du_dt.*conj(u_current) + conj(du_dt).*u_current;
end


% compute contribution to the derivative of the energy in each mode for
% each ROM at each resolution
for j = 1:length(N_list)
    
    % compute N and M for the current ROM
    N = N_list(j);
    M = 3*N;
    
    % construct output arrays of the proper size
    exact = zeros(2*N,2*N,2,t);
    R0 = zeros(2*N,2*N,2,t);
    R1 = zeros(2*N,2*N,2,t);
    R2 = zeros(2*N,2*N,2,t);
    R3 = zeros(2*N,2*N,2,t);
    R4 = zeros(2*N,2*N,2,t);
    
    % compute the indices corresponding to terms in this ROM from the exact
    % data computed above
    res_exact = [1:N,2*M_full - N+1:2*M_full];
    
    % make k array for the ROM
    k_vec = [0:M-1,-M:1:-1];
    [kx,ky] = ndgrid(k_vec,k_vec);
    k = zeros(2*M,2*M,2);
    k(:,:,1) = kx;
    k(:,:,2) = ky;
    
    % compute indices of resolved and unresolved modes
    a = 2:M;
    b = 2*M:-1:M+2;
    a_tilde = N+1:M;
    a_tilde2 = 2*N+1:M;
    
    % loop through every timestep
    for i = 1:t
        fprintf('Current resolution is N = %i\n',N)
        current_t = t_list(i);
        fprintf('Current time is t = %i out of %i\n\n',current_t,t_list(end))
        % find the exact derivative for the terms in this ROM
        exact_full = u_fullify(squeeze(exact_everything(:,:,:,:,i)),M_full);
        exact(:,:,:,i) = exact_full(res_exact,res_exact,:);
        
        % extract u and convert it into a full array with the proper
        % elements eliminated
        temp_u = squeeze(u(:,:,:,:,i));
        temp_u_full = u_fullify(temp_u,M_full);
        u_current = u_squishify(temp_u_full,N);
        u_full = u_fullify(u_current,M);
        
        % compute the R_k^i terms for i = 1,2,3,4
        [~,t0hat,t0tilde] = markov_term(u_full,a,b,k,a_tilde,a_tilde2);
        t0 = u_squishify(t0hat,N);
        % compute the energy derivative due to the R_0 term
        t0_energy = t0.*conj(u_current) + conj(t0).*u_current;
        t0_energy = u_fullify(t0_energy,M_full);
        R0(:,:,:,i) = t0_energy(res_exact,res_exact,:);
        
        [~,t1hat,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde,a_tilde2);
        t1 = u_squishify(t1hat,N);
        if time
            t1 = t1*current_t;
        end
        % compute the energy derivative due to the R_1 term
        t1_energy = t1.*conj(u_current) + conj(t1).*u_current;
        t1_energy = u_fullify(t1_energy,M_full);
         R1(:,:,:,i) = t1_energy(res_exact,res_exact,:);
         
         
         [t2,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2);
         t2 = u_squishify(t2,N);
         if time
             t2 = t2*current_t^2;
         end
         % compute the energy derivative due to the R_2 term
         t2_energy = t2.*conj(u_current) + conj(t2).*u_current;
         t2_energy = u_fullify(t2_energy,M_full);
         R2(:,:,:,i) = t2_energy(res_exact,res_exact,:);
         
         
         
         [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde,a_tilde2);
         t3 = u_squishify(t3,N);
         if time
             t3 = t3*current_t^3;
         end
         % compute the energy derivative due to the R_3 term
         t3_energy = t3.*conj(u_current) + conj(t3).*u_current;
         t3_energy = u_fullify(t3_energy,M_full);
         R3(:,:,:,i) = t3_energy(res_exact,res_exact,:);
         
         
         
         t4 = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde,a_tilde2);
         t4 = u_squishify(t4,N);
         if time
             t4 = t4*current_t^4;
         end
         % compute the energy derivative due to the R_4 term
         t4_energy = t4.*conj(u_current) + conj(t4).*u_current;
         t4_energy = u_fullify(t4_energy,M_full);
         R4(:,:,:,i) = t4_energy(res_exact,res_exact,:);
        
        
        
    end
    
    % compute the RHS for the least squares solve
    RHS = R0 - exact;
    err_list = squeeze(max(abs(RHS),[],[1,2,3]));
    err_list_exact = squeeze(max(abs(exact),[],[1,2,3]));
    err_list0 = squeeze(max(abs(R0),[],[1,2,3]));
    err_list1 = squeeze(max(abs(R1),[],[1,2,3]));
    err_list2 = squeeze(max(abs(R2),[],[1,2,3]));
    err_list3 = squeeze(max(abs(R3),[],[1,2,3]));
    err_list4 = squeeze(max(abs(R4),[],[1,2,3]));
    figure(1)
    xlabel("Time")
    ylabel("Max norm of diff. between Markov term and exact sol.")
    title("2D Euler exact - Markov error")
    hold on
    plot(t_list, err_list)
    figure(2)
    title("R1")
    xlabel("Time")
    ylabel("Max norm of R1")
    hold on
    plot(t_list, err_list1)
    figure(3)
    title("R2")
    xlabel("Time")
    ylabel("Max norm of R2")
    hold on
    plot(t_list, err_list2)
    figure(4)
    title("R3")
    xlabel("Time")
    ylabel("Max norm of R3")
    hold on
    plot(t_list, err_list3)
    figure(5)
    title("R4")
    xlabel("Time")
    ylabel("Max norm of R4")
    hold on
    plot(t_list, err_list4)
    figure(6)
    title("R0")
    xlabel("Time")
    ylabel("Max norm of R0")
    hold on
    plot(t_list,err_list0)
    figure(7)
    title("exact solution")
    xlabel("Time")
    ylabel("Max norm of exact solution")
    hold on
    plot(t_list, err_list_exact)
    Legend{j} = sprintf('N=%d',N) 
    drawnow

end
figure(1)
legend(Legend,'Location','northwest','FontSize',10)
figure(2)
legend(Legend,'Location','northwest','FontSize',10)
figure(3)
legend(Legend,'Location','northwest','FontSize',10)
figure(4)
legend(Legend,'Location','northwest','FontSize',10)
figure(5)
legend(Legend,'Location','northwest','FontSize',10)
figure(6)
legend(Legend,'Location','northwest','FontSize',10)
figure(7)
legend(Legend,'Location','northwest','FontSize',10)
