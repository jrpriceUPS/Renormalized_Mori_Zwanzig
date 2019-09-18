function v = decay_init(N)
%
%

u = rand(2*N,2*N,2);
u = fftn_norm(u);

% create k
kvec2 = [0:N-1,-N:1:-1];
[kx,ky] = ndgrid(kvec2,kvec2);
k = zeros(2*N,2*N,2);
k(:,:,1) = kx;
k(:,:,2) = ky;

v = zeros(2*N,2*N,2);

for i = 1:2*N
    for j = 1:2*N
        if i ~= 1 && j ~= 1
            decay = exp(norm(squeeze(k(i,j,:))));
            v(i,j,:) = u(i,j,:)./decay;
        end
    end
end

v = incomp_init(v);

