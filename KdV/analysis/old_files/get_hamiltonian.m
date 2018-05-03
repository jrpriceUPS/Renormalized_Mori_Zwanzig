function ham = get_hamiltonian(u,N,epsilon)

a = size(u);

ux = 1j*u.*repmat((0:a(1)-1).',1,a(2));


u_full = zeros(a(1)*3,a(2));
ux_full = zeros(a(1)*3,a(2));

u_full(1:N,:) = u(1:N,:);
u_full(a(1)*3-N+2:a(1)*3,:) = flipud(conj(u(2:N,:)));

ux_full(1:N,:) = ux(1:N,:);
ux_full(a(1)*3-N+2:a(1)*3,:) = flipud(conj(ux(2:N,:)));

u_real = ifft_norm(u_full);
ux_real = ifft_norm(ux_full);


ham = sum(u_real.^3 - 3 * epsilon^2*ux_real.^2,1);