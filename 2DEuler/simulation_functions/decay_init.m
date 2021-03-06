function v = decay_init(u_full)
%
%
M = length(u_full);
N = M/2;
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
            v(i,j,:) = u_full(i,j,:)./decay^2./norm(squeeze(u_full(i,j,:)));
        end
    end
end


