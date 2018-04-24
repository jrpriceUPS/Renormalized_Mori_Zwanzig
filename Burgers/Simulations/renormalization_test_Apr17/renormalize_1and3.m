function [c1,c3] = renormalize_1and3(alpha,N_list,u_list,t_list,exact_derivative,time)

markov_energy = zeros(max(N_list),length(t_list));
tmodel_energy = zeros(max(N_list),length(t_list));
t3model_energy = zeros(max(N_list),length(t_list));

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
        
        %compute t^3-model term
        [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde);
        t3model_energy(1:N,i,j) = t3(1:N).*conj(u) + conj(t3(1:N)).*u;
        if time
            t3model_energy(1:N,i,j) = t3model_energy(1:N,i,j)*t^3;
        end
    end
    
end


c1 = zeros(1,length(N_list));
c3 = zeros(2,length(N_list));

for i = 1:length(N_list)
%     if time
%         window = abs(sum(squeeze(tmodel_energy(1:N,i,:)),1)) > 1e-16 & abs(sum(squeeze(tmodel_energy(1:N,i,:)),1)) < 1e-10;
%     else
%         window = t_list.*abs(sum(squeeze(tmodel_energy(1:N,i,:)),1)) > 1e-16 & t_list.*abs(sum(squeeze(tmodel_energy(1:N,i,:)),1)) < 1e-10;
%     end
    window = 1:length(t_list);
    
    N = N_list(i);
    exact = exact_derivative(1:N,window);
    R0 = squeeze(markov_energy(1:N,i,window));
    R1 = squeeze(tmodel_energy(1:N,i,window));
    R3 = squeeze(t3model_energy(1:N,i,window));
    
    % compute the RHS for the least squares solve
    RHS = R0 - exact;
    
    b = -sum(RHS(:).*R1(:));
    
    % construct the matrix for the least squares solve
    A11 = sum(R1(:).*R1(:));
    
    c1(1,i) = A11\b;
    
    
    %         figure(1)
    %         subplot(2,2,1)
    %         hold off
    %         plot(t_list(window),sum(exact,1),'b')
    %         hold on
    %         plot(t_list(window),sum(R0 + R1*c1(1,i,j),1),'r')
    
    
    
    
    
    
    
    
    b = [b
        -sum(RHS(:).*R3(:))];
    
    % construct the matrix for the least squares solve
    A13 = sum(R1(:).*R3(:));
    
    A31 = A13;
    A33 = sum(R3(:).*R3(:));
    
    % solve the system and store the result
    A = [A11 A13
         A31 A33];
    
    c3(1:2,i) = A\b;
    
    %         figure(1)
    %         subplot(2,2,3)
    %         hold off
    %         plot(t_list(window),sum(exact,1),'b')
    %         hold on
    %         plot(t_list(window),sum(R0 + R1*c3(1,i,j) + R3*c3(2,i,j),1),'r')
    
    
    
end

figure
subplot(2,1,1)
plot(log(N_list),log(c1),'b.','markersize',20)
hold on
plot(log(N_list),log(c3(1,:)),'k.','markersize',20)
xlabel('log(N)','fontsize',16)
ylabel('log(a)','fontsize',16)
title('t-model coefficient','fontsize',16)
legend('Only fit t-model','fit t-model and t^3-model')

subplot(2,1,2)
plot(log(N_list),log(-c3(2,:)),'k.','markersize',20)
title('t^3-model coefficient','fontsize',16)
legend('fit t-model and t^3-model')
xlabel('log(N)','fontsize',16)
ylabel('log(-a)','fontsize',16)



