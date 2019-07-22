function u = taylor_green(X,Y)
%
%The Taylor-Green vortex initial condition
%
% u = [ sin(x)cos(y)
%      -cos(x)sin(y)]
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
% X  =  NxN array of x indices as produced by ndgrid
%
% Y  =  NxN array of y indices as produced by ndgrid
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% u  =  an NxNx2 array of evaluations of the Taylor-Green vortex
%       u(i,j,a) gives the x_i, y_j, entry for dimension a


u = zeros(length(X),length(X),2);
u(:,:,1) = sin(X).*cos(Y);
u(:,:,2) = -cos(X).*sin(Y);