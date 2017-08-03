function [u_rec,t_rec]=fullsolution_ground(M,epsilon,endtime,dt,rec,tol)
%
%[u_rec,t_rec]=fullsolution_ground(M,epsilon,endtime,dt,rec,tol)
%
%Computes the "full" solution of KdV equation with a pseudo-fourth order
%implicit-explicit Runge-Kutta solver. Uses timestep dt, but only records
%data every rec steps.
%
%u_t + 6*u*u_x + epsilon^2*u_{xxx} = 0
%Periodic boundary conditions, initial condition u(x,0)=sin(x)
%
%Solved in Fourier space with modes [-M,M].
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
%Example Input: [u_rec,t_rec]=fullsolution_ground(64,0.1,1,1e-4,1000,1e-10);
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
B=@(x,t) nonlinear_full(x);

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