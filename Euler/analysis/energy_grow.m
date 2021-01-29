function [value,isterminal,direction] = energy_grow(~,u,init_energy,percent_increase,params)

N = params.N;
M = params.M;

u = reshape(u,[N,N,N,3,4]);

% isolate modes that do not include any components larger than N
a = [1:N,2*M-N+1:2*M];
    
u_full = u_fullify(u,M);
current_energy = sum(sum(sum(sum(u_full(a,a,a,:).*conj(u_full(a,a,a,:))))));

value = double(current_energy < init_energy*(1+percent_increase));
if isnan(current_energy)
   value = 0;
end
isterminal = 1;
direction = 0;
end