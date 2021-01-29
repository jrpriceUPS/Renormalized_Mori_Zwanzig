function energy = get_3D_energy(u,N)
%
%Calculates the energy in the first all modes where every entry of the
%wavevector is N or smaller
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  M x M x M x 3 x 4 x length(t_list) array of Fourier modes at all timesteps
%
%  N  =  maximal mode to include when computing energy
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  energy  =  the energy in modes 1:N at all times

s = size(u);
M = s(1);

% isolate modes that do not include any components larger than N
a = [1:N,2*M-N+1:2*M];

t = s(end);
energy = zeros(t,1);

for i = 1:t
    
    u_full = u_fullify(u(:,:,:,:,:,i),M);
    energy(i) = (1/2)*sum(sum(sum(sum(u_full(a,a,a,:).*conj(u_full(a,a,a,:))))));

end