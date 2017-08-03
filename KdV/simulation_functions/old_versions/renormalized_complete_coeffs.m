function [t_list,u_list]=renormalized_complete_coeffs(N,epsilon,alpha,endtime,dt,coeffs,howoften,tol,blowup)
%
%[t_list,u_list]=renormalized_complete(N,epsilon,alpha,endtime,dt,howoften,tol1,tol2,blowup)
%
%Computes the complete renormalized ROM solution for N modes based on a 
%"full" model of 2N modes of KdV equation with a pseudo-fourth order 
%implicit-explicit Runge-Kutta solver.
%
%u_t + alpha*u*u_x + epsilon^2*u_{xxx} = 0
%Periodic boundary conditions, initial condition u(x,0)=sin(x)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       N  =  number of positive Fourier modes in t-model
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
%    tol1  =  tolerance for violation of conservation of energy
%
%    tol2  =  tolerance for when to compute renormalization coefficients
%
%  blowup  =  logical variable indicating how to handle loss of energy
%             conservation (if 1, save the data up to current time, if 0 
%             give error)
%
%Example Input: [t_list,u_list]=renormalized_complete(16,0.1,1,10,1e-3,1,1e-10,1e-10,0)
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
x=linspace(0,2*pi*(4*N-1)/(4*N),4*N);

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

%identify the number of steps that need to be taken
nsteps = length(t)-1;

F_modes = [1:N,2*N:4*N+2,5*N+2:6*N];
G_modes = N+1:5*N+1;
k = [0:3*N-1,-3*N:-1].';

%compute cutoff for different regimes of behavior
cutoff = ceil((2.8/(dt*epsilon^2))^(1/3));

if cutoff>N
    ns  = 1:N;
    s = [];
else
    ns  = 1:cutoff;
    s = cutoff+1:N;
end

%define the linear and nonlinear portions of the right hand side
A=diag(1j*((1:N)-1).^3*epsilon^2);
B=@(x,t) renormalized_complete_3rd_order(x,3*N,N,epsilon,alpha,F_modes,G_modes,k,t,coeffs);

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%construct vector of modes we will advance (just the positive modes), and
%fill it
u = zeros(N,1);
u(:) = u_complete(1:N);

%save the initial total energy
total_energy = sum(2*abs(u).^2);

%create u storage matrix
u_list = zeros(N,length(t_list));
u_list(:,1) = u(1:N);
current_index = 1;


%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

for i=1:nsteps
    
    u = RK4_stiff_nonstiff_step(u,t(i),dt,A,B,s,ns,As);
    
    %if total energy in increases above the initial amount, abort the evaluation
    if (sum(2*abs(u).^2)-total_energy)/total_energy>tol
        
        if blowup
            u_list = u_list(:,1:current_index);
            t_list = t_list(1:current_index);
            break
        else
            
            error('Energy decay of t-model not preserved, try a smaller timestep')
            
        end
    end
    
    %save solution periodically
    if mod(i,howoften)==0
        current_index = current_index + 1;
        u_list(:,current_index) = u;
    end
end