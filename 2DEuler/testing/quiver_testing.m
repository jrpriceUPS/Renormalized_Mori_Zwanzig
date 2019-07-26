
clear all;close all;
N = 12;
M = 3*N;

% uniform grid
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N).';
y = x;
    
 % 2D array of data points
[X,Y] = ndgrid(x,y);

eval = taylor_green(X,Y);
%u_full = fftn_norm(eval);
size(eval);

x = X;
y = Y;

%[x,y] = meshgrid(linspace(0,2*pi*(2*M-1)/(2*M),2*M),linspace(0,2*pi*(2*M-1)/(2*M),2*M));
%renorm_u = u_full./sqrt((u_full(:,:,1).^2 + u_full(:,:,2)^2));
u = eval(:,:,1);
v = eval(:,:,2);

%quiverC2D(x,y,u,v,1000,.5);
plot_vfield(u,v,N,1000,.5);

% [x,y] = meshgrid(linspace(0,10,100),linspace(0,10,100));
%        u = exp(-0.2*(x-5).^2 - 0.2*(y-5).^2);
%        v = -u;
%        quiverC2D(x,y,u,v,1000);