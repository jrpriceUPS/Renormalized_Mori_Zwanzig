function v = incomp_init(u_full)
% 
% Random initial condition that satisfies incompressibility
%
% div(u) = 0
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%       u_full  =  a 2Nx2Nx2 array
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% v  =  a 2Nx2Nx2 array values that satisfies div(v)=0
%       v(i,j,a) gives the x_i, y_j, entry for dimension a

% get the size
M = length(u_full);
N = M/2;

% create k
kvec2 = [0:N-1,-N:1:-1];
[kx,ky] = ndgrid(kvec2,kvec2);
k = zeros(2*N,2*N,2);
k(:,:,1) = kx;
k(:,:,2) = ky;

% fill out v = (I - (kk')/(|k|^2))u_k
v = zeros(M,M,2);
for i = 1:M
    for j = 1:M
        if i ~= 1 && j ~= 1
            k_vec = squeeze(k(i,j,:));
            u_vec = squeeze(u_full(i,j,:));
            v(i,j,:) = (eye(2)-((k_vec*k_vec')/norm(k_vec)^2))*u_vec;
        end        
    end
end




