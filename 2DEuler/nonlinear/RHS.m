function du_dt = RHS(u,time,params)
%
% The ODE of an RHS of an ROM of 2D Euler's equations:
%
% du_k/dt = func(u)
%
%
%%%%%%%%%
%INPUTS:%
%%%%%%%%%
%
%      u  =  the current state array (NxNx2 unspooled into a vector for ode45)
%
%   time  =  current time in simulation
%
%
% The structure params contains the other necessary parameters for
% simulation
%
%      k  =  the array of wavevectors (NxNx2)
%
%      N  =  the resolution of the model
%
%      M  =  the resolution of the full model (ROM) / array size for dealiasing (full)
%
%   func  =  the equation for the RHS (for example, t4model_RHS)
%
%  coeff  =  constant coefficient assigned to t-model
%
%
%%%%%%%%%%
%OUTPUTS:%
%%%%%%%%%%
%
% du_dt  =  the time derivative of the state vector (NxNx2 unspooled into
%           a vector for ode45)

if params.print_time
   time 
end

N = params.N;
M = params.M;

u = reshape(u,[N,N,2,2]);

% find the full state of u
params.u_full = u_fullify(u,M);

params.time = time;

% fill the appropriate entries of du_dt with the full euler RHS
du_dt = params.func(params);

% unspool the derivative so it can be used by ode45
du_dt = du_dt(:);

