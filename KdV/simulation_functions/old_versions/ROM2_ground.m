function [u_rec,t_rec]=ROM2_ground(M,epsilon,endtime,dt,rec,tol)
%
%[u_rec,t_rec]=tmodel_ground(M,epsilon,endtime,dt,tol)
%
%Computes the 2nd order ROM solution for M modes based on a "full" model of 2M
%modes of KdV equation with a pseudo-fourth order implicit-explicit 
%Runge-Kutta solver. Uses timestep dt, but only records data every rec 
%steps.
%
%u_t + 6*u*u_x + epsilon^2*u_{xxx} = 0
%Periodic boundary conditions, initial condition u(x,0)=sin(x)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       M  =  number of positive Fourier modes
%
% epsilon  =  degree of dispersion
%
% endtime  =  time at which simulation ends
%
%      dt  =  timestep
%
%     rec  =  amount of data points to record (evenly distributed)
%
%     tol  =  tolerance for violation of conservation of energy
%
%
%Example Input: [u_rec,t_rec]=ROM2_ground(64,0.1,1,1e-4,1000,1e-10);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_rec  =  a matrix where each column corresponds to the full fft of the
%            solution at the timestep of the associated entry in t_rec
%
%  t_rec  =  a vector of the recording times in the simulation

%%%%%%%%%%%%%%%%%%%
%Initial Condition%
%%%%%%%%%%%%%%%%%%%

%define the ordinates in real space with sufficient buffer
x=linspace(0,2*pi*(2*M-1)/(2*M),2*M);

%define the initial condition as the Fourier transform of the sine function
%and ensure that high energy modes are zero in spite of any numerical 
%errors
u_complete=fft(sin(x)).';

%save the initial total energy
total_energy = sum(abs(u_complete/(2*M)).^2);

%identify indexes of different Fourier modes of interest
positivemodes = 1:M;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Time-Stepping Variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct vector of times
t=0:dt:endtime;

%construct vector of times of recording
t_rec = t(1:1/(dt*rec):length(t));

%construct output matrix for data
u_rec = zeros(length(x),length(t_rec));
u_rec(:,1) = u_complete/(2*M);

%define the linear and nonlinear portions of the right hand side
A=diag(1j*(positivemodes-1).^3*epsilon^2);
B=@(x,t) nonlinear_2nd_order(x,2*M,epsilon,t);

%identify the number of steps that need to be taken
nsteps = length(t)-1;

%compute cutoff for different regimes of behavior
cutoff = ceil((2.8/(dt*epsilon^2))^(1/3));

if cutoff>M   
    ns  = 1:M;
    s = [];   
else 
    ns  = 1:cutoff;
    s = cutoff+1:M;
end

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%construct vector of modes we will advance (just the positive modes), and
%fill it
u = zeros(M,1);
u(:) = u_complete(positivemodes);

%initialize record indicator
record = 2;


%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

for i=1:nsteps
%     %compute U_1
%     U(:,1) = u;
%     
%     %compute linear and nonlinear parts of RHS from U_1
%     nonlin(:,1) = B(U(:,1));
%     linear(:,1) = A*U(:,1);
%     
%     %compute U_2
%     U(ns,2) = u(ns)+1/2*dt*(linear(ns,1)+nonlin(ns,1));   
%     if ~isempty(s)
%        U(s,2) = (eye(length(s))-1/3*dt*As)\(u(s)+1/2*dt*nonlin(s,1)+1/6*dt*linear(s,1)); 
%     end
%     
%     %compute linear and nonlinear parts of RHS from U_2
%     nonlin(:,2) = B(U(:,2));
%     linear(:,2) = A*U(:,2);
%     
%     %compute U_3
%     U(ns,3) = u(ns)+1/2*dt*(linear(ns,2)+nonlin(ns,2));
%     if ~isempty(s)
%         U(s,3) = (eye(length(s))-dt*As)\(u(s)+1/2*dt*nonlin(s,2)+1/2*dt*linear(s,1)-dt*linear(s,2));
%     end
%     
%     %compute linear and nonlinear parts of RHS from U_3
%     nonlin(:,3) = B(U(:,3));
%     linear(:,3) = A*U(:,3);
%     
%     %compute U_4
%     U(ns,4) = u(ns)+dt*(linear(ns,3)+nonlin(ns,3));
%     if ~isempty(s)
%         U(s,4) = (eye(length(s))-1/3*dt*As)\(u(s)+dt*nonlin(s,3)+2/3*dt*linear(s,3));
%     end
%     
%     %compute linear and nonlinear parts of RHS from U_4
%     nonlin(:,4) = B(U(:,4));
%     linear(:,4) = A*U(:,4);
%     
%     %use U_1, U_2, U_3, and U_4 to update complete a Runge-Kutta timestep
%     u(ns) = u(ns) + 1/6*dt*(linear(ns,1)+nonlin(ns,1)+linear(ns,4)+nonlin(ns,4))...
%                   + 1/3*dt*(linear(ns,2)+nonlin(ns,2)+linear(ns,3)+nonlin(ns,3));
%               
%     if ~isempty(s)
%         u(s) = u(s) + 1/6*dt*(linear(s,1)+nonlin(s,1)+linear(s,4)+nonlin(s,4))...
%                     + 1/3*dt*(linear(s,2)+nonlin(s,2)+linear(s,3)+nonlin(s,3));
%     end

    u = RK4_stiff_nonstiff_step(u,t(i),dt,A,B,s,ns,As);

    if t(i+1)==t_rec(record)
        u_rec(1:M,record) = u/(2*M);
        u_rec(M+2:2*M,record) = conj(flipud(u(2:M)))/(2*M);
        record=record+1;
    end
    
    %if total energy conservation fails, print an error and quit
    if (sum(abs(u_rec(:,record-1)/(2*M)).^2)-total_energy)/total_energy>tol
        error('Energy conservation violated. Try a smaller timestep.')
    end
end