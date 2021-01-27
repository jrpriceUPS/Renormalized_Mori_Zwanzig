function [R0,R1,R2,R3,R4,RTS,RT] = ROM_dE_memory(c4,N,t_list,u_list,alpha)
%
% [R0,R1,R2,R3,R4,RT] = ROM_dE_memory(c4,N,t_list,u_list,alpha)
%
% Computes the contribution to the change in energy for each renormalized memory term
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%            c4  =  5 x 1 array of renormalization coefficients and tau for
%            specified N
%
%            N  =  size of the ROM solution
%
%            u_list  =  N x length(t_list) array of the ROM solution
%
%            t_list  =  the times associated with the u_list solutions
%
%            alpha  =  degree of dispersion
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
%       RT = Total change in energy of all orders
%
%       RTS = Integral of the total change in energy of all orders,
%       integrated over the entire time domain (primarily a check
%       condition)

% Preallocate arrays
markov_energy = zeros(N,length(t_list));
tmodel_energy = zeros(N,length(t_list));
t2model_energy = zeros(N,length(t_list));
t3model_energy = zeros(N,length(t_list));
t4model_energy = zeros(N,length(t_list));

R0 = zeros(1,length(t_list));
R1 = zeros(1,length(t_list));
R2 = zeros(1,length(t_list));
R3 = zeros(1,length(t_list));
R4 = zeros(1,length(t_list));
RT = zeros(1,length(t_list));

F_modes = [1:N,2*N+1:4*N+1,5*N+1:6*N];
G_modes = [N+1:5*N+1];
M = 3*N;

for i = 1:length(t_list)
    
    disp(sprintf('Calculating energy for N = %i, time t = %i',N,t_list(i)))
    t = t_list(i);
    u = u_list(1:N,i);
    
    % Markov term
    [t0,t0hat,t0tilde,u_full] = markov_term_Burgers(u,M,N,alpha,F_modes,G_modes);
    markov_energy(1:N,i) = t0(1:N).*conj(u) + conj(t0(1:N)).*u;
    R0(1,i) = sum(markov_energy(1:N,i),1);
    
    % t-model term
    [t1,~,~] = tmodel_term_Burgers(u_full,t0tilde,alpha,F_modes,G_modes);
    tmodel_energy(1:N,i) = c4(1)*t^(1-1*c4(5))*(t1(1:N).*conj(u) + conj(t1(1:N)).*u);
    R1(1,i) = sum(tmodel_energy(1:N,i),1);
    
    % t^2-model term
    [t2,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde] = t2model_term_Burgers(u_full,alpha,t0hat,t0tilde,F_modes,G_modes);
    t2model_energy(1:N,i) = c4(2)*t^(2-2*c4(5))*(t2(1:N).*conj(u) + conj(t2(1:N)).*u);
    R2(1,i) = sum(t2model_energy(1:N,i),1);
    
    % t^3-model term
    [t3,Ehat,Etilde,Fhat,Ftilde] = t3model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dtilde);
    t3model_energy(1:N,i) = c4(3)*t^(3-3*c4(5))*(t3(1:N).*conj(u) + conj(t3(1:N)).*u);
    R3(1,i) = sum(t3model_energy(1:N,i),1);
    
    % t^4-model term
    t4 = t4model_term_Burgers(alpha,F_modes,G_modes,u_full,t0hat,t0tilde,Ahat,Atilde,Bhat,Btilde,Dhat,Dtilde,Ehat,Etilde,Fhat,Ftilde);
    t4model_energy(1:N,i) = c4(4)*t^(4-4*c4(5))*(t4(1:N).*conj(u) + conj(t4(1:N)).*u);
    R4(1,i) = sum(t4model_energy(1:N,i),1);

    % Total energy change
    RT(1,i) = R0(1,i)+R1(1,i)+R2(1,i)+R3(1,i)+R4(1,i);

end

RTS = trapz(t_list,RT.');