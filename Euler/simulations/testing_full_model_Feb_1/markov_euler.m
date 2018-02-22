function du_dt = markov_euler(u,k,N)
%
% The ODE for the markov term of 3D Euler's equations:
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
% after projection, this becomes
%
% du_k/dt = Ck(u_hat,u_hat)
%
% where u_hat is a vector with only the resolved modes
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

% full model has twice as many modes in each direction
M = 2*N;

u = reshape(u,[N,N,N,3,4]);

% compute indices of modes
a = 2:N;
b = 2*M:-1:2*M-N+2;

% resolved modes
a_hat = 2:N;
a_tilde = N+1:M;

% find the convolution of uu.'
u_full = u_fullify(u,M);

% fill the appropriate entries of du_dt with the full euler RHS
du_dt = markov_term(u_full,a,b,k,a_hat,a_tilde);

% unspool the derivative so it can be used by ode45
du_dt = du_dt(:);

