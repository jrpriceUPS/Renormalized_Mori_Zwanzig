function [t_list,u_list]=fullsolution_4(M,epsilon,alpha,endtime,dt,howoften,tol)
%
%Input: [t_list,u_list]=fullsolution_4(M,epsilon,alpha,endtime,dt,howoften,tol)
%
%Computes the "full" solution of KdV equation with a pseudo-fourth order
%implicit-explicit Runge-Kutta solver
%
%u_t + alpha*u*u_x + epsilon^2*u_{xxx} = 0
%Periodic boundary conditions, initial condition u(x,0)=sin(x)
%
%Solved in Fourier space with modes [-M,M].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       M  =  number of positive Fourier modes in t-model
%
% epsilon  =  degree of dispersion
%
%   alpha  =  coefficient of nonlinearity
%
% endtime  =  time at which simulation ends
%
%      dt  =  timestep
%
%howoften  =  how frequently to save the solution (i.e. 10 -> every 10
%             timesteps)
%
%     tol  =  tolerance for violation of conservation of energy
%
%
%Example Input: [t_list,u_list]=fullsolution_4(64,0.1,6,1,1e-4,100,1e-10)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t_list  =  vector of the times at which data were recorded
%
%  u_list  =  an Mxlength(t_list) matrix containing the soltion at each
%             time in t_list

%%%%%%%%%%%%%%%%%%%
%Initial Condition%
%%%%%%%%%%%%%%%%%%%

%define the ordinates in real space with sufficient buffer
x=linspace(0,2*pi*(2*M-1)/(2*M),2*M);

%define the initial condition as the Fourier transform of the sine function
%and ensure that high energy modes are zero in spite of any numerical
%errors
u_complete=fft_norm(sin(x).');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Time-Stepping Variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct vector of times
t=0:dt:endtime;
t_list = 0:dt*howoften:endtime;

%define the linear and nonlinear portions of the right hand side
A=diag(1j*((1:M)-1).^3*epsilon^2);
B=@(x,t) nonlinear_full(x,M,alpha);

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
u(:) = u_complete(1:M);

%compute total energy
total_energy = sum(2*abs(u).^2);

%create u storage matrix
u_list = zeros(length(u),length(t_list));
u_list(:,1) = u;
current_index = 1;


%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

for i=1:nsteps
    u = RK4_stiff_nonstiff_step_old(u,t(i),dt,A,B,s,ns,As);

    %if total energy conservation fails, print an error and quit
    if abs(sum(2*abs(u).^2)-total_energy)/abs(total_energy)>tol
        error('Energy conservation violated. Try a smaller timestep.')
    end
    
    %plot the real space solution periodically
    if mod(i,howoften)==0
        current_index = current_index + 1;
        u_list(:,current_index) = u;
    end
end