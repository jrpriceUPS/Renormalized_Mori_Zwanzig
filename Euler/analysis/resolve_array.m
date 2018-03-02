function [u_trim,t_trim,tmodel_size_list] = resolve_array(u,t,tol)

s = size(u);
N = s(1);
M = 3*N;

a = 2:M;
b = 2*M:-1:M+2;
a_tilde = N+1:M;

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

t_index = s(6);

t_model_size = 0;
i = 0;

tmodel_size_list = zeros(t_index,1);

while t_model_size < tol
    i = i + 1;
    current_u_temp = squeeze(u(:,:,:,:,:,i));
    
    u_full = u_fullify(current_u_temp,M);
    time = t(i)
    
    [~,~,t0tilde] = markov_term(u_full,a,b,k,a_tilde);
    [~,t1hat,~] = tmodel_term(u_full,t0tilde,a,b,k,a_tilde);
    
    t_model_size = sum(abs(t1hat(:)).^2)
    tmodel_size_list(i) = t_model_size;
end

u_trim = u(:,:,:,:,:,1:i);
t_trim = t(1:i);