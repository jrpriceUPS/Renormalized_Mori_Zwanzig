function du_dt = full_euler(u,k,N,t)
%
% The ODE for the full 3D Euler's equations:
%
% u_t + u dot grad(u) = -grad(p), grad dot u = 0
%
% In Fourier space:
%
% du_k/dt = -i sum_{p+q=k} k dot u_p A_k u_q
% where A_k = I - kk.'/|k|^2
%
% The convolutions are handled by transforming into real space. So the real
% ODE becomes:
%
% du_k/dt = -i [FFT(uu.')]_k k
%           + ik/|k|^2(k.'[FFT(uu.')]_k k)
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
% u  =  the current state array (NxNxNx3 unspooled into a vector for ode45)
%
% k  =  the array of wavevectors (NxNxNx3)
%
% N  =  the resolution of the model
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% du_dt  =  the time derivative of the state vector (NxNxNx3 unspooled into
%           a vector for ode45)

M = 3*N;

u = reshape(u,[N,N,N,3,4]);

%compute indices of modes
a = 2:N;
b = 2*M:-1:2*M-N+2;

% find the convolution of uu.'
u_full = u_fullify(u,M);

% fill the appropriate entries of du_dt with the full euler RHS
du_dt = full_RHS(u_full,a,b,k,N);

% unspool the derivative so it can be used by ode45
du_dt = du_dt(:);

