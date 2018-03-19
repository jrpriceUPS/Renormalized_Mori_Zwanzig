%clear all;close all;

% load relevant folders into the path
addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

N = 8;
M = 3*N;

% uniform grid
x = linspace(0,2*pi*(2*M-1)/(2*M),2*M).';
y = x;
z = x;

% 3D array of data points
[X,Y,Z] = ndgrid(x,y,z);

% create initial condition
eval = taylor_green(X,Y,Z);
u_full = fftn_norm(eval);

% make k array
k_vec = [0:M-1,-M:1:-1];
[kx,ky,kz] = ndgrid(k_vec,k_vec,k_vec);
k = zeros(2*M,2*M,2*M,3);
k(:,:,:,1) = kx;
k(:,:,:,2) = ky;
k(:,:,:,3) = kz;

a = 2:M;
b = 2*M:-1:M+2;
a_tilde = N+1:M;

tic
C1 = Ck(u_full,u_full,a,b,k,a_tilde);
toc

tic
C2 = Ck_old(u_full,u_full,a,b,k,a_tilde);
toc

max(abs(C1(:)-C2(:)))