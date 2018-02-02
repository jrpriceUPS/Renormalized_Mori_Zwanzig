function u = taylor_green(X,Y,Z)
%
%The Taylor-Green vortex initial condition
%
% u = [ sin(x)cos(y)cos(z)
%      -cos(x)sin(y)cos(z)
%      0]
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
% X  =  NxNxN array of x indices as produced by ndgrid
%
% Y  =  NxNxN array of y indices as produced by ndgrid
%
% Z  =  NxNxN array of z indices as produced by ndgrid
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% u  =  an NxNxNx3 array of evaluations of the Taylor-Green vortex
%       u(i,j,k,a) gives the x_i, y_j, z_k entry for dimension a


u = zeros(length(X),length(X),length(X),3);
u(:,:,:,1) = sin(X).*cos(Y).*cos(Z);
u(:,:,:,2) = -cos(X).*sin(Y).*cos(Z);
u(:,:,:,3) = 0;