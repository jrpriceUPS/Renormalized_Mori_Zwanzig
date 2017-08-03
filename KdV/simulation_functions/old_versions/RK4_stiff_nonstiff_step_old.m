function u = RK4_stiff_nonstiff_step_old(u,t,dt,A,B,s,ns,As)
%
%u = RK4_stiff_nonstiff_step(u,dt,A,B,s,ns)
%
%Computes a single RK4 step using an implicit/explicit scheme to handle
%stiff and non-stiff parts of the problem
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       u  =  the positive Fourier modes ([0:N]) that we seek to advance
%
%       t  =  current time
%
%      dt  =  timestep
%
%       A  =  matrix for linear part of problem
%
%       B  =  nonlinear part of problem in form of @(x,t) B(x,t)
%
%       s  =  indices of stiff modes
%
%      ns  =  indices of non-stiff modes
%
%      As  =  stiff part of the linear part of the problem
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u  =  u at the next timestep


%compute the length of u
M = length(u);

%initialize data structure containing nonlinear and linear right hand side
%values for intermediate steps
nonlin = zeros(M,4);
linear = zeros(M,4);

%initialize data structure containing intermediate steps
U = zeros(M,4);


%compute U_1
U(:,1) = u;

%compute linear and nonlinear parts of RHS from U_1
nonlin(:,1) = B(U(:,1),t);
linear(:,1) = A*U(:,1);

%compute U_2
U(ns,2) = u(ns)+1/2*dt*(linear(ns,1)+nonlin(ns,1));
if ~isempty(s)
    U(s,2) = (eye(length(s))-1/3*dt*As)\(u(s)+1/2*dt*nonlin(s,1)+1/6*dt*linear(s,1));
end

%compute linear and nonlinear parts of RHS from U_2
nonlin(:,2) = B(U(:,2),t+dt/2);
linear(:,2) = A*U(:,2);

%compute U_3
U(ns,3) = u(ns)+1/2*dt*(linear(ns,2)+nonlin(ns,2));
if ~isempty(s)
    U(s,3) = (eye(length(s))-dt*As)\(u(s)+1/2*dt*nonlin(s,2)+1/2*dt*linear(s,1)-dt*linear(s,2));
end

%compute linear and nonlinear parts of RHS from U_3
nonlin(:,3) = B(U(:,3),t+dt/2);
linear(:,3) = A*U(:,3);

%compute U_4
U(ns,4) = u(ns)+dt*(linear(ns,3)+nonlin(ns,3));
if ~isempty(s)
    U(s,4) = (eye(length(s))-1/3*dt*As)\(u(s)+dt*nonlin(s,3)+2/3*dt*linear(s,3));
end

%compute linear and nonlinear parts of RHS from U_4
nonlin(:,4) = B(U(:,4),t+dt);
linear(:,4) = A*U(:,4);

%use U_1, U_2, U_3, and U_4 to update complete a Runge-Kutta timestep
u(ns) = u(ns) + 1/6*dt*(linear(ns,1)+nonlin(ns,1)+linear(ns,4)+nonlin(ns,4))...
    + 1/3*dt*(linear(ns,2)+nonlin(ns,2)+linear(ns,3)+nonlin(ns,3));

if ~isempty(s)
    u(s) = u(s) + 1/6*dt*(linear(s,1)+nonlin(s,1)+linear(s,4)+nonlin(s,4))...
        + 1/3*dt*(linear(s,2)+nonlin(s,2)+linear(s,3)+nonlin(s,3));
end