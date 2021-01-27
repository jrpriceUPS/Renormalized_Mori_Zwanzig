function [R0,R1,R2,R3,R4,RTS,RT] = ROM_dE_memory(c4,N,t_list,u)
    %
    % [R0,R1,R2,R3,R4,RTS,RT] = ROM_dE_memory(N,t_list,u_list)
    %
    % Computes the contribution to the change in energy for non-renormalized t-model term
    %
    %
    %%%%%%%%%
    %INPUTS:%
    %%%%%%%%%
    %
    %           c4  =  5 x 1 array of renormalization coefficients and tau for
    %           specified N
    %
    %            N  =  size of the ROM solution
    %
    %            u_list  =  NxNxNx3x4xlength(t_list) array of the ROM solution
    %
    %            t_list  =  the times associated with the u_list solutions
    %
    %
    %%%%%%%%%%
    %OUTPUTS:%
    %%%%%%%%%%
    % 
    %       R0 = 1 x length(t_list) for change in energy of the Markov term summed over all modes -> should be zero!
    %
    %       R1 = 1 x length(t_list)for change in energy of the t-model summed over all modes
    %
    %       R2 = 1 x length(t_list)for change in energy of the t2-model summed over all modes
    %
    %       R3 = 1 x length(t_list)for change in energy of the t3-model summed over all modes
    %
    %       R4 = 1 x length(t_list)for change in energy of the t4-model summed over all modes
    %
    %       RT = total change in energy of all orders
    %
    %       RTS = Integral of the total change in energy of all orders,
    %       integrated over the entire time domain (primarily a check
    %       condition)

    % Preallocate arrays 
    R0 = zeros(1,length(t_list));
    R1 = zeros(1,length(t_list));
    R2 = zeros(1,length(t_list));
    R3 = zeros(1,length(t_list));
    R4 = zeros(1,length(t_list));
    RT = zeros(1,length(t_list));
    
    M = 3*N;

    % Make k array for the ROM
    k_vec = [0:M-1,-M:1:-1];
    [kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
    k = zeros(2*M,2*M,2*M,3);
    k(:,:,:,1) = kx;
    k(:,:,:,2) = ky;
    k(:,:,:,3) = kz;
    
    % Compute indices of resolved and unresolved modes
    a = 2:M;
    b = 2*M:-1:M+2;
    a_tilde = N+1:M;
    a_tilde2 = 2*N+1:M;

    % isolate modes that do not include any components larger than N
    res = [1:N,2*M-N+1:2*M];

    for i = 1:length(t_list)

        current_t = t_list(i);

        disp(sprintf('Current time is t = %i out of %i\n',current_t,t_list(end)))

        % extract u and convert it into a full array with the proper
        % elements eliminated
        u_current = squeeze(u(:,:,:,:,:,i));
     
        u_full = u_fullify(u_current,M);
 
        [~,t0hat,t0tilde] = markov_term(u_full,a,b,k,a_tilde,a_tilde2);
        t0 = u_squishify(t0hat,N);

        % compute the energy derivative due to the R_0 term
        t0_energy = t0.*conj(u_current) + conj(t0).*u_current;
        t0_energy_full = u_fullify(t0_energy,M);
        R0(1,i) = (1/2)*sum(t0_energy_full(res,res,res,:),[1 2 3 4]);

        [~,t1hat,t1tilde] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde,a_tilde2);
        t1 = u_squishify(t1hat,N);
        
        % compute the energy derivative due to the R_1 term
        t1_energy = c4(1)*current_t^(1-1*c4(5))*(t1.*conj(u_current) + conj(t1).*u_current);
        t1_energy_full = u_fullify(t1_energy,M);
        R1(1,i) = (1/2)*sum(t1_energy_full(res,res,res,:),[1 2 3 4]);

        [t2hat,Ahat,Atilde,Bhat,Btilde] = t2model_term(u_full,t0hat,t0tilde,t1tilde,a,b,k,a_tilde,a_tilde2);
        t2 = u_squishify(t2hat,N);
        
        % compute the energy derivative due to the R_2 term
        t2_energy = c4(2)*current_t^(2-2*c4(5))*(t2.*conj(u_current) + conj(t2).*u_current);
        t2_energy_full = u_fullify(t2_energy,M);
        R2(1,i) = (1/2)*sum(t2_energy_full(res,res,res,:),[1 2 3 4]);

        [t3hat,Ehat,Etilde,Fhat,Ftilde] = t3model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Btilde,a,b,k,a_tilde,a_tilde2);
        t3 = u_squishify(t3hat,N);
        
        % compute the energy derivative due to the R_3 term
        t3_energy = c4(3)*current_t^(3-3*c4(5))*(t3.*conj(u_current) + conj(t3).*u_current);
        t3_energy_full = u_fullify(t3_energy,M);
        R3(1,i) = (1/2)*sum(t3_energy_full(res,res,res,:),[1 2 3 4]);

        t4hat = t4model_term(u_full,t0hat,t0tilde,t1hat,t1tilde,Ahat,Atilde,Bhat,Btilde,Ehat,Etilde,Fhat,Ftilde,a,b,k,a_tilde,a_tilde2);
        t4 = u_squishify(t4hat,N);
        
        % compute the energy derivative due to the R_4 term
        t4_energy = c4(4)*current_t^(4-4*c4(5))*(t4.*conj(u_current) + conj(t4).*u_current);
        t4_energy_full = u_fullify(t4_energy,M);
        R4(1,i) = (1/2)*sum(t4_energy_full(res,res,res,:),[1 2 3 4]);

        RT(1,i) = (R0(1,i)+R1(1,i)+R2(1,i)+R3(1,i)+R4(1,i)).';

    end

    RTS = trapz(t_list,RT);

