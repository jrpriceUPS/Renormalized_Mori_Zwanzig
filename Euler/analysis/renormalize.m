function coeff_array = renormalize(u,N_list,t_list,time)

s = size(u);
N_full = s(1);
M_full = 3*N_full;
t = s(6);

exact_everything = zeros(s);

coeff_array = zeros(4,length(N_list));


% make k array
k_vec_full = [0:M_full-1,-M_full:1:-1];
[kx_full,ky_full,kz_full] = ndgrid(k_vec_full,k_vec_full,k_vec_full);
k_full = zeros(2*M_full,2*M_full,2*M_full,3);
k_full(:,:,:,1) = kx_full;
k_full(:,:,:,2) = ky_full;
k_full(:,:,:,3) = kz_full;

a_full = 2:M_full;
b_full = 2*M_full:-1:M_full+2;

for i = 1:2
    
    u_current = squeeze(u(:,:,:,:,:,i));
    
    u_full = u_fullify(u_current,M_full);
    
    du_dt = Ck(u_full,u_full,a_full,b_full,k_full,[]);
    du_dt = u_squishify(du_dt,N_full);
    
    exact_everything(:,:,:,:,:,i) = du_dt.*conj(u_current) + conj(du_dt).*u_current;
    
end





for j = 1:length(N_list);
    N = N_list(j)
    M = 2*N;
    
    exact = zeros(M,M,M,3,t);
    R0 = zeros(M,M,M,3,t);
    R1 = zeros(M,M,M,3,t);
    R2 = zeros(M,M,M,3,t);
    R3 = zeros(M,M,M,3,t);
    R4 = zeros(M,M,M,3,t);
    

    res_exact = [1:N,2*M_full - N+1:2*M_full];
    
    % make k array
    k_vec = [0:M-1,-M:1:-1];
    [kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
    k = zeros(2*M,2*M,2*M,3);
    k(:,:,:,1) = kx;
    k(:,:,:,2) = ky;
    k(:,:,:,3) = kz;
    
    a = 2:M;
    b = 2*M:-1:M+2;
    a_tilde = N+1:M;
    
    for i = 1:t
        current_t = t_list(i)

        exact_full = u_fullify(squeeze(exact_everything(:,:,:,:,:,i)),M_full);
        exact(:,:,:,:,i) = exact_full(res_exact,res_exact,res_exact,:);
        
        temp_u = squeeze(u(:,:,:,:,:,i));
        temp_u_full = u_fullify(temp_u,M_full);
        u_current = u_squishify(temp_u_full,N);
        u_full = u_fullify(u,M);
        
        [~,t0hat,t0tilde] = markov_term(u_full,a,b,k,a_tilde);
        t0 = u_squishify(t0hat,N);
        [~,t1hat,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde);
        t1 = u_squishify(t1hat,N);
        [t2,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde);
        t2 = u_squishify(t2,N);
        [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde);
        t3 = u_squishify(t3,N);
        t4 = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde);
        t4 = u_squishify(t4,N);
        
        if time
            
            t1 = t1*current_t;
            t2 = t2*current_t^2;
            t3 = t3*current_t^3;
            t4 = t4*current_t^4;
            
        end
        
        
        t0_energy = t0.*conj(u_current) + conj(t0).*u_current;
        t0_energy = u_fullify(t0_energy,M_full);
        R0(:,:,:,:,i) = t0_energy(res_exact,res_exact,res_exact,:);
        
        
        t1_energy = t1.*conj(u_current) + conj(t1).*u_current;
        t1_energy = u_fullify(t1_energy,M_full);
        R1(:,:,:,:,i) = t1_energy(res_exact,res_exact,res_exact,:);
        
        
        t2_energy = t2.*conj(u_current) + conj(t2).*u_current;
        t2_energy = u_fullify(t2_energy,M_full);
        R2(:,:,:,:,i) = t2_energy(res_exact,res_exact,res_exact,:);
        
        
        t3_energy = t3.*conj(u_current) + conj(t3).*u_current;
        t3_energy = u_fullify(t3_energy,M_full);
        R3(:,:,:,:,i) = t3_energy(res_exact,res_exact,res_exact,:);
        
        
        t4_energy = t4.*conj(u_current) + conj(t4).*u_current;
        t4_energy = u_fullify(t4_energy,M_full);
        R4(:,:,:,:,i) = t4_energy(res_exact,res_exact,res_exact,:);
        
    end
    
    R0 = R0 - exact;
    
    b = -[sum(R0(:).*R1(:))
          sum(R0(:).*R2(:))
          sum(R0(:).*R3(:))
          sum(R0(:).*R4(:))];
     
     A11 = sum(R1(:).*R1(:));
     A12 = sum(R1(:).*R2(:));
     A13 = sum(R1(:).*R3(:));
     A14 = sum(R1(:).*R4(:));
     
     A21 = A12;
     A22 = sum(R2(:).*R2(:));
     A23 = sum(R2(:).*R3(:));
     A24 = sum(R2(:).*R4(:));
     
     A31 = A13;
     A32 = A23;
     A33 = sum(R3(:).*R3(:));
     A34 = sum(R3(:).*R4(:));
     
     A41 = A14;
     A42 = A24;
     A43 = A34;
     A44 = sum(R4(:).*R4(:));
    
    
     A = [A11 A12 A13 A14
          A21 A22 A23 A24
          A31 A32 A33 A34
          A41 A42 A43 A44];
      
    coeff_array(:,j) = A\b
end
