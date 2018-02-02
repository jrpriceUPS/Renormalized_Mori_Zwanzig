% old code used to test the effectiveness of u_fullify and u_squishify

%clear all;close all;
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N = 4;

%test_func = rand(2*N,2*N,2*N,3)-0.5;
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

%test_func = taylor_green(X,Y,Z);
test_func = rand([size(X),3]);

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