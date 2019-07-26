N = 4;
M = 8;

k_vec = [0:M-1,-M:-1];
[kx,ky] = ndgrid(k_vec,k_vec);
k = zeros(2*M,2*M,2);
k(:,:,1) = kx;
k(:,:,2) = ky;

a = 2:M;
b = 2*M:-1:M+2;
a_tilde = N+1:M;
a_tilde2 = 2*N+1:M;

u = rand(2*M,2*M,2);
v = rand(2*M,2*M,2);

u_full = fftn_norm2(u);
v_full = fftn_norm2(v);

test1 = Ck(u_full,v_full,a,b,k,a_tilde,a_tilde2);
test2 = Ck(u_full,u_full,a,b,k,a_tilde,a_tilde2);
test3 = Ck(u_full,5*u_full + 3*v_full,a,b,k,a_tilde,a_tilde2);

max(abs(test3(:) - 3*test1(:)-5*test2(:)));