function [u_trim,tmodel_size_list] = resolve_array(u,tol)

s = size(u);
N = s(1);
M = 3*N;

a = 2:M_full;
b = 2*M_full:-1:M_full+2;
a_tilde = N_full+1:M_full;

t_index = s(6);

t_model_size = 0;
i = 0;

tmodel_size_list = zeros(t_index,1);

while t_model_size < tol
    
    i = i + 1;
    current_u_temp = squeeze(u(:,:,:,:,:,i));
    
    u_full = u_fullify(current_u_temp,M);
    
    [~,~,t0tilde] = markov_term(u_full,a,b,k,a_tilde);
    [~,t1hat,~] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde);
    
    t_model_size = sum(abs(t1hat(:)).^2);
    tmodel_size_list(i) = t_model_size;
    
end

u_trim = u(:,:,:,:,:,1:i);