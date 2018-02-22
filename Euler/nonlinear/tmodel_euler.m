function du_dt = tmodel_euler(u,k,N,time,coeff)
%
% The ODE for the t-model term of 3D Euler's equations:
%
% du_k/dt = Dk(u_hat,Ctilde(u_hat,u_hat))
%
% where u_hat is a vector with only the resolved modes
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%      u  =  the current state array (NxNxNx3 unspooled into a vector for ode45)
%
%      k  =  the array of wavevectors (NxNxNx3)
%
%      N  =  the resolution of the model
%
%   time  =  current time in simulation
%
%  coeff  =  constant coefficient assigned to t-model
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
a = 2:M;
b = 2*M:-1:M+2;

% unresolved modes
a_tilde = N+1:M;

% find the full state of u
u_full = u_fullify(u,M);

% fill the appropriate entries of du_dt with the full euler RHS
du_dt = tmodel_RHS(u_full,a,b,k,a_tilde,N,time,coeff);

% unspool the derivative so it can be used by ode45
du_dt = du_dt(:);

