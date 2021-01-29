function energy = get_energy(u,N)
%
%Calculates the energy in the first N modes at each time step for the 
%solution matrix u
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  N x length(t_list) array of Fourier modes at all timesteps
%
%  N  =  maximal mode to include when computing energy
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  energy  =  the energy following definition of 1/2 sum_{k \in F} |u_k|^2

s = size(u);
t = s(2);

u_full = zeros(2*N,t);
u_full(1:N,:) = u(1:N,:);
u_full(N+2:2*N,:) = conj(flipud(u(2:N,:)));

energy = (1/2)*sum(u_full.*conj(u_full)).';