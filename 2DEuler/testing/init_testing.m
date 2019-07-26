N = 2;

eval_input = fftn_norm(10*rand(2*N,2*N,2));
%eval_input = u_fullify(u_squishify(eval_input,N),M);
eval = incomp_init(eval_input);

% create k
kvec2 = [0:N-1,-N:1:-1];
[kx,ky] = ndgrid(kvec2,kvec2);
k = zeros(2*N,2*N,2);
k(:,:,1) = kx;
k(:,:,2) = ky;

sum = 0;
for i = 1:2*N
    for j = 1:2*N
        k_vec = squeeze(k(i,j,:));
        eval_vec = squeeze(eval(i,j,:));
        sum = sum + dot(k_vec,eval_vec);
    end
end

sum

%sum
M = 2*N;


% create k
kvec2 = [0:M-1,-M:1:-1];
[kx,ky] = ndgrid(kvec2,kvec2);
k = zeros(2*M,2*M,2);
k(:,:,1) = kx;
k(:,:,2) = ky;


% load relevant parameters into parameter structure
params.k = k;
params.N = N;
params.M = M;
params.func = @(x) full_RHS(x);
params.coeff = [];
params.a = 2:M;
params.b = 2*M:-1:M+2;
params.a_tilde = N+1:M;
params.a_tilde2 = 2*N+1:M;
params.print_time = 1;

%eval = u_squishify(eval,N);
% 
% v_deriv = RHS(eval(:),0,params);
% v_deriv_new = zeros(N,N,2,2,2);
% v_deriv_new(:,:,:,:,1) = reshape(v_deriv,[N,N,2,2]);
% get_2D_energy(v_deriv_new,N);

u_fullify(u_squishify(eval,N),M)
eval
