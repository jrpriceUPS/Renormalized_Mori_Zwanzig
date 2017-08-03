function [t_list,u_list,coeffs]=renormalized_4old2(N,epsilon,alpha,endtime,dt,howoften,tol,tol2)
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

%use different resolutions to compute data for renormalization coefficient
%fitting
N_list = [N,N-1,N-2,N-3];

nonlinear_terms = zeros(4,3);
nonlinear_markov = zeros(4,1);
time_recorded = zeros(4,1);

for i = 1:4
    
    N0 = N_list(i);
    
    %compute information about modes
    F_modes = [1:N0,2*N0+1:2*N0+3,3*N0+4:4*N0+2];
    G_modes = N0+1:3*N0+3;
    k = [0:2*N0,-2*N0-1:-1].';
    
    A_red = diag(1j*((1:N0)-1).^3*epsilon^2);
    B_red = @(x,t) renormalized_3rd_order(x,2*N0+1,N0,epsilon,alpha,F_modes,G_modes,k,t,[1;1;1]);
    
    if cutoff>N0
        ns_red  = 1:N0;
        s_red = [];
    else
        ns_red  = 1:cutoff;
        s_red = cutoff+1:N0;
    end
    
    %isolate the stiff portion of the linear portion of the right hand side
    As_red = A_red(s_red,s_red);
    
    u_red = u(1:N0);
    
    
    threshold = 0;
    t = 0;
    num_steps = 0;
    
    
    while threshold < tol2
        t = t + dt;
        num_steps = num_steps + 1;
        
        u_red = RK4_stiff_nonstiff_step(u_red,t,dt,A_red,B_red,s_red,ns_red,As_red);
        
        %compute the total L2 norm being drained from resolved modes
        [nonlin0,u_full] = markov_term(u_red,2*N0+1,N0,alpha);
        
        %compute t-model term
        [nonlin1,uu_star] = tmodel_term(u_full,nonlin0,alpha,F_modes);
        threshold = abs(t*sum(nonlin1(1:N0).*conj(u_red) + conj(nonlin1(1:N0)).*u_red));
        
    end
    
    %compute t^2-model term
    [nonlin2,uk,uu,uk_uu_u,uk_uu_u_star] = t2model_term(u_full,nonlin0,uu_star,alpha,F_modes,G_modes,k,epsilon);
    
    %compute t^3-model term
    nonlin3 = t3model_term(alpha,F_modes,k,epsilon,u_full,uu,uu_star,uk,uk_uu_u,uk_uu_u_star);
    
    nonlinear_markov(i) = sum(nonlin0(1:N0).*conj(u_red) + conj(nonlin0(1:N0)).*u_red);
    nonlinear_terms(i,1) = t*sum(nonlin1(1:N0).*conj(u_red) + conj(nonlin1(1:N0)).*u_red);
    nonlinear_terms(i,2) = -t^2/2*sum(nonlin2(1:N0).*conj(u_red) + conj(nonlin2(1:N0)).*u_red);
    nonlinear_terms(i,3) = t^3/6*sum(nonlin3(1:N0).*conj(u_red) + conj(nonlin3(1:N0 )).*u_red);
    
    time_recorded(i) = num_steps;
    
end


%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

jettison = 0;
how_many_done = 0;
RHS_match = zeros(4,1);

for i=1:nsteps
    
    u = RK4_stiff_nonstiff_step(u,t_list(i),dt,A,B,s,ns,As);
    
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
        else
            u_list(:,current_index) = u;
        end
    end
    
    if jettison == 0
        
        time_to_match = find(i+1 == time_recorded);
        
        if ~isempty(time_to_match)
            
            how_many_done = how_many_done + length(time_to_match);
            [RHS,~] = markov_term(u,3*N,2*N,alpha);
            
            for j = 1:length(time_to_match)
                
                which_ROM = time_to_match(j);
                RHS_match(which_ROM) = sum(RHS(1:N_list(which_ROM)).*conj(u(1:N_list(which_ROM)))+conj(RHS(1:N_list(which_ROM))).*u(1:N_list(which_ROM)));
                RHS_match(which_ROM) = RHS_match(which_ROM) - nonlinear_markov(which_ROM);
                
            end
            
            if how_many_done == 4
                jettison = 1;
                u = u(1:N);
                coeffs = coeffs_calculator_new(RHS_match,nonlinear_terms,N_list.')
                
                F_modes = [1:N,2*N+1:2*N+3,3*N+4:4*N+2];
                G_modes = N+1:3*N+3;
                k = [0:2*N,-2*N-1:-1].';
                A = A(1:N,1:N);
                if cutoff>N
                    ns  = 1:N;
                    s = [];
                else
                    ns  = 1:cutoff;
                    s = cutoff+1:N;
                end
                As = A(s,s);
                B = @(x,t) renormalized_3rd_order(x,2*N+1,N,epsilon,alpha,F_modes,G_modes,k,t,coeffs);
            end
            
        end
    end
    
    
    
end

coeffs = coeffs.';