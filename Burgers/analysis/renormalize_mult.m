function [c1,c2,c3,c4] = renormalize_mult(alpha,N_list,u_list,t_list,exact_derivative,time)

markov_energy = zeros(max(N_list),length(t_list));
tmodel_energy = zeros(max(N_list),length(t_list));
t2model_energy = zeros(max(N_list),length(t_list));
t3model_energy = zeros(max(N_list),length(t_list));
t4model_energy = zeros(max(N_list),length(t_list));

for i = 1:length(N_list)
    
    N = N_list(i);
    %gather parameters needed for simulation
    F_modes = [1:N,2*N:4*N+2,5*N+2:6*N];
    G_modes = N+1:5*N+1;
    M = 3*N;
    
    for j = 1:length(t_list)
        
        disp(sprintf('Currently calculating coefficients for N = %i, time t = %i',N,t_list(j)))
        t = t_list(j);
        u = u_list(1:N,j);
        
        %compute Markov term
        [t0,t0hat,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
        markov_energy(1:N,i,j) = t0(1:N).*conj(u) + conj(t0(1:N)).*u;
        
        %compute t-model term
        [t1,~,~] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes);
        tmodel_energy(1:N,i,j) = t1(1:N).*conj(u) + conj(t1(1:N)).*u;
        if time
            tmodel_energy(1:N,i,j) = tmodel_energy(1:N,i,j)*t;
        end
        
        %compute t^2-model term
        [t2,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde] = t2model_term_Burgers(u_full,alpha,t0hat,t0tilde,F_modes,G_modes);
        t2model_energy(1:N,i,j) = t2(1:N).*conj(u) + conj(t2(1:N)).*u;
        if time
            t2model_energy(1:N,i,j) = t2model_energy(1:N,i,j)*t^2;
        end
        
        %compute t^3-model term
        [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde);
        t3model_energy(1:N,i,j) = t3(1:N).*conj(u) + conj(t3(1:N)).*u;
        if time
            t3model_energy(1:N,i,j) = t3model_energy(1:N,i,j)*t^3;
        end
        
        
        %compute t^4-model term
        t4 = t4model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde,Ehat,Etilde,Fhat,Ftilde);
        t4model_energy(1:N,i,j) = t4(1:N).*conj(u) + conj(t4(1:N)).*u;
        if time
            t4model_energy(1:N,i,j) = t4model_energy(1:N,i,j)*t^4;
        end
    end
    
end


c1 = zeros(1,length(N_list),length(t_list));
c2 = zeros(2,length(N_list),length(t_list));
c3 = zeros(3,length(N_list),length(t_list));
c4 = zeros(4,length(N_list),length(t_list));

for j = 1:length(t_list)
    for i = 1:length(N_list)
        window = 1:j;
        
        N = N_list(i);
        exact = exact_derivative(1:N,window);
        R0 = squeeze(markov_energy(1:N,i,window));
        R1 = squeeze(tmodel_energy(1:N,i,window));
        R2 = squeeze(t2model_energy(1:N,i,window));
        R3 = squeeze(t3model_energy(1:N,i,window));
        R4 = squeeze(t4model_energy(1:N,i,window));
        
        % compute the RHS for the least squares solve
        RHS = R0 - exact;
        
        b = -sum(RHS(:).*R1(:));
        
        % construct the matrix for the least squares solve
        A11 = sum(R1(:).*R1(:));
        
        c1(1,i,j) = A11\b;
        
        
        %         figure(1)
        %         subplot(2,2,1)
        %         hold off
        %         plot(t_list(window),sum(exact,1),'b')
        %         hold on
        %         plot(t_list(window),sum(R0 + R1*c1(1,i,j),1),'r')
        
        
        
        
        b = [b
            -sum(RHS(:).*R2(:))];
        
        % construct the matrix for the least squares solve
        A12 = sum(R1(:).*R2(:));
        
        A21 = A12;
        A22 = sum(R2(:).*R2(:));
        
        % solve the system and store the result
        A = [A11 A12
            A21 A22];
        
        c2(1:2,i,j) = A\b;
        
        %         figure(1)
        %         subplot(2,2,2)
        %         hold off
        %         plot(t_list(window),sum(exact,1),'b')
        %         hold on
        %         plot(t_list(window),sum(R0 + R1*c2(1,i,j) + R2*c2(2,i,j),1),'r')
        
        
        
        
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
        
        c3(1:3,i,j) = A\b;
        
        %         figure(1)
        %         subplot(2,2,3)
        %         hold off
        %         plot(t_list(window),sum(exact,1),'b')
        %         hold on
        %         plot(t_list(window),sum(R0 + R1*c3(1,i,j) + R2*c3(2,i,j) + R3*c3(3,i,j),1),'r')
        
        
        
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
        
        c4(1:4,i,j) = A\b;
        
        %         figure(1)
        %         subplot(2,2,4)
        %         hold off
        %         plot(t_list(window),sum(exact,1),'b')
        %         hold on
        %         plot(t_list(window),sum(R0 + R1*c4(1,i,j) + R2*c4(2,i,j) + R3*c4(3,i,j) + R4*c4(4,i,j),1),'r')
        
        
    end
    
%     if mod(j,10) == 0
%         figure(1)
%         subplot(2,2,1)
%         hold off
%         plot(log(N_list),log(c1(1,:,j)),'b.','markersize',20)
%         hold on
%         plot(log(N_list),log(c2(1,:,j)),'r.','markersize',20)
%         plot(log(N_list),log(c3(1,:,j)),'k.','markersize',20)
%         plot(log(N_list),log(c4(1,:,j)),'c.','markersize',20)
%         xlabel('log(N)','fontsize',16)
%         ylabel('log(a)','fontsize',16)
%         title(sprintf('t-model coefficient at t = %i',t_list(j)),'fontsize',16)
%         legend('fit ROM 1','fit ROM 1+2','fit ROM 1+2+3','fit ROM 1+2+3+4')
%         
%         subplot(2,2,2)
%         hold off
%         plot(log(N_list),log(c2(2,:,j)),'r.','markersize',20)
%         hold on
%         plot(log(N_list),log(c3(2,:,j)),'k.','markersize',20)
%         plot(log(N_list),log(c4(2,:,j)),'c.','markersize',20)
%         xlabel('log(N)','fontsize',16)
%         ylabel('log(a)','fontsize',16)
%         title(sprintf('t^2-model coefficient at t = %i',t_list(j)),'fontsize',16)
%         legend('fit ROM 1+2','fit ROM 1+2+3','fit ROM 1+2+3+4')
%         
%         subplot(2,2,3)
%         hold off
%         plot(log(N_list),log(c3(3,:,j)),'k.','markersize',20)
%         hold on
%         plot(log(N_list),log(c4(3,:,j)),'c.','markersize',20)
%         xlabel('log(N)','fontsize',16)
%         ylabel('log(a)','fontsize',16)
%         title(sprintf('t^3-model coefficient at t = %i',t_list(j)),'fontsize',16)
%         legend('fit ROM 1+2+3','fit ROM 1+2+3+4')
%         
%         subplot(2,2,4)
%         hold off
%         plot(log(N_list),log(c4(4,:,j)),'c.','markersize',20)
%         xlabel('log(N)','fontsize',16)
%         ylabel('log(a)','fontsize',16)
%         title(sprintf('t^4-model coefficient at t = %i',t_list(j)),'fontsize',16)
%         legend('fit ROM 1+2+3+4')
%     end
    
end




