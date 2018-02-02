% an old algorithm used to test how u_squishify and u_fullify worked

%clear all;close all;

N = 4;

%test_func = rand(2*N,2*N,2*N,3)-0.5;
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

test_func = taylor_green(X,Y,Z);

u_full = fftn_norm(test_func);
u_full(N+1,:,:,:) = 0;
u_full(:,N+1,:,:) = 0;
u_full(:,:,N+1,:) = 0;


u = u_squishify(u_full,N);
u_full2 = u_fullify(u,2*N);
u2 = u_squishify(u_full2,N);



x = u2;
y = u;

max(max(max(max(max(abs(x-y))))))