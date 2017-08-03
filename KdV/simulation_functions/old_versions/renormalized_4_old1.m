function [t_list,u_list,coeffs,resid,u_end,t_norm]=renormalized_4_old1(N,epsilon,alpha,endtime,dt,howoften,tol,tol2)
%
%[t_list,u_list]=renormalized_4(N,epsilon,endtime,dt,howoften,tol)
%
%Computes the t-model solution for N modes based on a "full" model of 2N
%modes of KdV equation with a pseudo-fourth order implicit-explicit
%Runge-Kutta solver.
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
%     tol  =  tolerance for violation of conservation of energy
%
%
%Example Input: [t_list,u_list]=tmodel_4(64,0.1,6,1,1e-4,100,1e-10)
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

%compute cutoff for different regimes of behavior
cutoff = ceil((2.8/(dt*epsilon^2))^(1/3));

if cutoff>2*N
    ns  = 1:2*N;
    s = [];
else
    ns  = 1:cutoff;
    s = cutoff+1:2*N;
end

%define the linear and nonlinear portions of the right hand side
A=diag(1j*((1:2*N)-1).^3*epsilon^2);
B=@(x,t) nonlinear_full(x,2*N,alpha);

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%construct vector of modes we will advance (just the positive modes), and
%fill it
u = zeros(2*N,1);
u(:) = u_complete(1:2*N);

%save the initial total energy
total_energy = sum(2*abs(u).^2);

%create u storage matrix
u_list = zeros(N,length(t_list));
u_list(:,1) = u(1:N);
current_index = 1;

%create coefficient storage matrix
coeffs = zeros(3,length(t_list)-1);
resid = zeros(1,length(t_list)-1);
u_end = zeros(1,length(t_list)-1);
t_norm = zeros(1,length(t_list)-1);

jettison = 0;
current_tnorm = tol2/10;

F_modes = cell(4,1);
G_modes = cell(4,1);
k = cell(4,1);
A_red = cell(4,1);
B_red = cell(4,1);
ns_red = cell(4,1);
s_red = cell(4,1);
As_red = cell(4,1);
u_red = cell(4,1);


N_list = [N,N-2,N-4,N-6];

for i = 1:4
    
    N0 = N_list(i);
    
    %compute information about modes
    F_modes{i} = [1:N0,2*N0+1:2*N0+3,3*N0+4:4*N0+2];
    G_modes{i} = N0+1:3*N0+3;
    k{i} = [0:2*N0,-2*N0-1:-1].';
    
    A_red{i} = diag(1j*((1:N0)-1).^3*epsilon^2);
    B_red{i} = @(x,t) renormalized_3rd_order(x,2*N0+1,N0,epsilon,alpha,F_modes{i},G_modes{i},k{i},t,[1;1;1]);
    
    if cutoff>N0
        ns_red{i}  = 1:N0;
        s_red{i} = [];
    else
        ns_red{i}  = 1:cutoff;
        s_red{i} = cutoff+1:N0;
    end
    
    %isolate the stiff portion of the linear portion of the right hand side
    As_red{i} = A_red{i}(s_red{i},s_red{i});
    
    u_red{i} = u(1:N0);
    
end


%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

for i=1:nsteps
    
    u = RK4_stiff_nonstiff_step(u,t(i),dt,A,B,s,ns,As);
    
    if jettison == 0
        for j = 1:4
            u_red{j} = RK4_stiff_nonstiff_step(u_red{j},t(i),dt,A_red{j},B_red{j},s_red{j},ns_red{j},As_red{j});
        end
    end
    
    %if total energy in increases above the initial amount, abort the evaluation
    if (sum(2*abs(u).^2)-total_energy)/total_energy>tol
        %error('Energy decay of t-model not preserved, try a smaller timestep')
        u_list = u_list(:,1:current_index);
        t_list = t_list(1:current_index);
        break
    end
    
    %save coefficients periodically
    if mod(i,howoften)==0
        current_index = current_index + 1;
        
        if jettison == 0
            u_list(:,current_index) = u(1:N);
            [coeffs(:,current_index-1),resid(current_index-1),t_norm(current_index-1)] = coeffs_calculator(u,u_red,2*N+1,N,epsilon,alpha,F_modes,G_modes,k,t(i+1));
            u_end(current_index-1) = 2*abs(u(2*N)).^2;
            current_tnorm = t_norm(current_index-1);
        else
            u_list(:,current_index) = u;
        end
    end
    
    if jettison == 0
        %if 2*abs(u(2*N)).^2 > tol2
        if abs(current_tnorm)>tol2
            u = u(1:N);
            coeffs_final = coeffs(:,current_index-1)
            B = @(x,t) renormalized_3rd_order(x,2*N+1,N,epsilon,alpha,F_modes,G_modes,k,t,coeffs_final);
            jettison = 1;
            
            coeffs = coeffs(:,1:current_index-1);
            resid = resid(1:current_index-1);
            u_end = u_end(1:current_index-1);
        end
    end
    
end

coeffs = coeffs.';