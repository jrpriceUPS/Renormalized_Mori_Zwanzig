function [u_complete,t,energyvector]=tmodel_2(M,N,epsilon,endtime,dt,fig,howoften,tol)
%
%Example input: [u_complete,t,energyvector]=tmodel_2(M,N,epsilon,endtime,dt,fig,howoften)
%
%Computes the t-model solution for N modes based on a "full" model of 2M
%modes of KdV equation with a second order implicit-explicit solver. 
%Records the energy in the first N positive Fourier modes.
%
%u_t + 6*u*u_x + epsilon^2*u_{xxx} = 0
%Periodic boundary conditions, initial condition u(x,0)=sin(x)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       M  =  number of positive Fourier modes in t-model
%
%       N  =  number of positive Fourier modes in energy recording
%
% epsilon  =  degree of dispersion
%
% endtime  =  time at which simulation ends
%
%      dt  =  timestep
%
%     fig  =  boolean set to one if figure outputs are desired, zero
%             otherwise
%
%howoften  =  how frequently to plot the solution (i.e. 10 -> every 10
%             timesteps)
%
%     tol  =  tolerance for violation of conservation of energy
%
%
%Example Input: [u,t,energy]=tmodel_2(32,8,0.1,1,1e-5,1,1000,1e-4);
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  u_complete  =  vector of magnitude of energy in Fourier modes ordered 
%                  like [0:M-1,-M:-1]
%
%           t  =  A vector of the times in the simulation
%
%energyvector  =  A vector of the energy in the specified [-N,N-1] modes
%                 at each time step

%%%%%%%%%%%%%%%%%%%
%Initial Condition%
%%%%%%%%%%%%%%%%%%%

%define the ordinates in real space with sufficient buffer
x=linspace(0,2*pi*(2*M-1)/(2*M),2*M);

%define the initial condition as the Fourier transform of the sine function
u_complete=fft(sin(x)).';

%identify indexes of different Fourier modes of interest
positivemodes = 1:M;
negativemodes = M+2:2*M;
Findex=[1:N,2*M-N+2:2*M];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Time-Stepping Variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct vector of times
t=0:dt:endtime;

%identify the number of steps that need to be taken
nsteps = length(t)-1;

%define the linear and nonlinear portions of the right hand side
A=diag(1j*(positivemodes-1).^3*epsilon^2);
B=@(x,t) nonlinear_red(x,2*M,t);

%initialize data structure containing nonlinear and linear right hand side
%values for past four steps
nonlin = zeros(M,5);
linear = zeros(M,5);

%construct vector of modes we will advance (just the positive modes), and
%fill it
u = zeros(M,1);
u(:) = u_complete(positivemodes);

%initialize total energy vector
energyvector=zeros(1,nsteps+1);

%fill in initial values for each
energyvector(1)=sum(abs(u_complete(Findex)/(2*M)).^2);

%save the initial total energy
total_energy = sum(abs(u_complete/(2*M)).^2);


%%%%%%%%%%%%%%%%%
%Initialize Plot%
%%%%%%%%%%%%%%%%%

%if we are plotting, initialize the figure window and plot the initial
%condition
if fig==1
    figure(1)
    plot(x,real(ifft(u_complete)));
    axis([0,2*pi,-3,3])
    title(sprintf('t = 0'))
end


%%%%%%%%%%%%%%%%%%%%
%Initial Four Steps%
%%%%%%%%%%%%%%%%%%%%

%The second order scheme requires data from the previous four steps. We
%compute these with our pseudo-fourth order method

%initialize data structure containing intermediate steps
U = zeros(M,4);

%initialize data structure containing nonlinear and linear right hand side
%values for intermediate steps
nonlin1 = zeros(M,4);
linear1 = zeros(M,4);

%compute cutoff for different regimes of behavior (ns -> not stiff,
%s -> stiff)
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

%complete four steps of pseudo-fourth order scheme
for i=1:4
    %compute U_1
    U(:,1) = u;
    
    %compute linear and nonlinear parts of RHS from U_1
    nonlin(:,1) = B(U(:,1),t(i));
    linear(:,1) = A*U(:,1);
    
    %save data that will be used to jump-start main simulation
    nonlin(:,i)=nonlin1(:,1);
    linear(:,i)=linear1(:,1);
    
    %compute U_2
    U(ns,2) = u(ns)+1/2*dt*(linear1(ns,1)+nonlin1(ns,1));   
    if ~isempty(s)
       U(s,2) = (eye(length(s))-1/3*dt*As)\(u(s)+1/2*dt*nonlin1(s,1)+1/6*dt*linear1(s,1)); 
    end
    
    %compute linear and nonlinear parts of RHS from U_2
    nonlin1(:,2) = B(U(:,2),t(i)+dt/2);
    linear1(:,2) = A*U(:,2);
    
    %compute U_3
    U(ns,3) = u(ns)+1/2*dt*(linear1(ns,2)+nonlin1(ns,2));
    if ~isempty(s)
        U(s,3) = (eye(length(s))-dt*As)\(u(s)+1/2*dt*nonlin1(s,2)+1/2*dt*linear1(s,1)-dt*linear1(s,2));
    end
    
    %compute linear and nonlinear parts of RHS from U_3
    nonlin1(:,3) = B(U(:,3),t(i)+dt/2);
    linear1(:,3) = A*U(:,3);
    
    %compute U_4
    U(ns,4) = u(ns)+dt*(linear1(ns,3)+nonlin1(ns,3));
    if ~isempty(s)
        U(s,4) = (eye(length(s))-dt*As)\(u(s)+dt*nonlin1(s,3)+2/3*dt*linear1(s,3));
    end
    
    %compute linear and nonlinear parts of RHS from U_4
    nonlin1(:,4) = B(U(:,4),t(i)+dt);
    linear1(:,4) = A*U(:,4);
    
    %use U_1, U_2, U_3, and U_4 to update complete a Runge-Kutta timestep
    u(ns) = u(ns) + 1/6*dt*(linear1(ns,1)+nonlin1(ns,1)+linear1(ns,4)+nonlin1(ns,4))...
        + 1/3*dt*(linear1(ns,2)+nonlin1(ns,2)+linear1(ns,3)+nonlin1(ns,3));
    
    if ~isempty(s)
        u(s) = u(s) + 1/6*dt*(linear1(s,1)+nonlin1(s,1)+linear1(s,4)+nonlin1(s,4))...
            + 1/3*dt*(linear1(s,2)+nonlin1(s,2)+linear1(s,3)+nonlin1(s,3));
    end
    
    %construct full vector of Fourier modes and compute the energy
    u_complete(positivemodes) = u;
    u_complete(negativemodes) = conj(flipud(u(2:M)));
    energyvector(i+1)=sum(abs(u_complete(Findex)/(2*M)).^2);
end



%%%%%%%%%%%%%%%%%
%Main Simulation%
%%%%%%%%%%%%%%%%%

%compute cutoffs for different regimes of behavior
cutoff_1 = ceil((0.43/(dt*epsilon^2))^(1/3));
cutoff_2 = ceil((1.36/(dt*epsilon^2))^(1/3));

if cutoff_1>M
    reg1=1:M;
    reg2=[];
    reg3=[];
elseif cutoff_2>M
    reg1=1:cutoff_1;
    reg2=cutoff_1+1:M;
    reg3=[];
else
    reg1=1:cutoff_1;
    reg2=cutoff_1+1:cutoff_2;
    reg3=cutoff_2+1:M;
end

%isolate the linear part of the right hand side for trouble
%regions
A2 = A(reg2,reg2);
A3 = A(reg3,reg3);

%complete remainder of simulation
for i=5:nsteps
    
    %note index in which next step will be stored
    index=mod(i,5)+1;
    
    %adjust indices where current step and past four steps are located
    i_n = mod(index-2,5)+1;
    i_nm1 = mod(index-3,5)+1;
    i_nm2 = mod(index-4,5)+1;
    i_nm3 = mod(index-5,5)+1;
    i_nm4 = mod(index-6,5)+1;
    
    %update with nonlinear and linear parts of current step
    nonlin(:,i_n) = B(u,t(i));
    linear(:,i_n) = A*u(:);
    
    %Adams-Bashforth for both parts of small k regime
    u(reg1) = u(reg1) + dt*...
        (55/24*(nonlin(reg1,i_n)+linear(reg1,i_n))...
        -59/24*(nonlin(reg1,i_nm1)+linear(reg1,i_nm1))...
        +37/24*(nonlin(reg1,i_nm2)+linear(reg1,i_nm2))...
        -3/8*(nonlin(reg1,i_nm3)+linear(reg1,i_nm3)));
    
    %Adams-Bashforth for nonlinear part of medium k regime, Adams-Moulton for
    %linear part
    if ~isempty(reg2)
        u(reg2) = (eye(length(reg2))-475/1440*dt*A2)\(u(reg2)+dt*(...
            55/24*nonlin(reg2,i_n)...
            -59/24*nonlin(reg2,i_nm1)...
            +37/24*nonlin(reg2,i_nm2)...
            -3/8*nonlin(reg2,i_nm3)...
            +1427/1440*linear(reg2,i_n)...
            -133/240*linear(reg2,i_nm1)...
            +241/720*linear(reg2,i_nm2)...
            -173/1440*linear(reg2,i_nm3)...
            +3/160*linear(reg2,i_nm4)));
    end
    
    %Adams-Bashforth for nonlinear part of large k regime, Adams-Moulton for
    %linear part
    if ~isempty(reg3)
        u(reg3) = (eye(length(reg3))-3*dt/4*A3)\(u(reg3) + dt*(...
            55/24*nonlin(reg3,i_n)...
            -59/24*nonlin(reg3,i_nm1)...
            +37/24*nonlin(reg3,i_nm2)...
            -3/8*nonlin(reg3,i_nm3)...
            +1/4*linear(reg3,i_nm1)));
    end
    
    %construct full vector of Fourier modes and compute the energy
    u_complete(positivemodes) = u(:);
    u_complete(negativemodes) = conj(flipud(u(2:M)));
    energyvector(i+1)=sum(abs(u_complete(Findex)/(2*M)).^2);
    
    %if total energy in increases above the initial amount, abort the evaluation
    if (sum(abs(u_complete/(2*M)).^2)-total_energy)/total_energy>tol
       error('Energy decay of t-model not preserved, try a smaller timestep') 
    end
    
    %plot the real space solution periodically
    if fig==1
        if mod(i,howoften)==0
            figure(1)
            plot(x,real(ifft(u_complete)));
            axis([0,2*pi,-3,3])
            title(sprintf('t = %g', t(i+1)))
            pause(0.1)
        end
    end
end

u_complete = u_complete/(2*M);