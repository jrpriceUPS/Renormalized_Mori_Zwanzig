function [simulation_params] = full_init_KdV(simulation_params)
%
%[simulation_params] = full_init_KdV(simulation_params)
%
%Takes initial model and simulation_params structures and initializes them
%for a standard full simulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
%
%       N                  =  number of positive resolved modes
%
%       epsilon            =  degree of dispersion
%
%       dt                 =  timestep
%
%       alpha              =  degree of nonlinearity
%
%       initial_condition  =  initial condition as a function of x
%
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
%
%       M        =  size of the full model
%
%       F_modes  =  a cell array of indices for resolved modes in the full
%                   model
%
%       G_modes  =  a cell array of indices for unresolved modes in the full
%                   model
%
%       k        =  a cell array of wave numbers in full model
%
%       u        =  state vector of positive modes of full model
%
%       A        =  matrix of linear part of KdV model
%
%       B        =  function handle of nonlinear part of RHS in full model
%
%       As       =  stiff portion of A to be handled implicitly
%
%       s        =  indices of stiff portion of A
%
%       ns       =  indices of nonstiff portions of A
%
%       cutoff   =  cutoff for stiff portions of A

%create shorthand for ROM system size
N = simulation_params.N;

%general simulation parameters
epsilon = simulation_params.epsilon;
dt      = simulation_params.dt;


%define the ordinates in real space
x=linspace(0,2*pi*(2*N-1)/(2*N),2*N);

%define the initial condition as the Fourier transform of the sine function
%and ensure that high energy modes are zero in spite of any numerical
%errors
u_complete=fft_norm(simulation_params.initial_condition(x).');

%initialize cells indicating index information, and populate them
simulation_params.F_modes = 0;
simulation_params.G_modes = 0;
simulation_params.k = 0;
simulation_params.M = 3*N/2;

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
B=@(x,t) nonlinear_full_KdV(x,simulation_params.N,simulation_params.alpha);

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%construct vector of modes we will advance in full model (just the positive modes), and
%fill it
u = zeros(N,1);
u(:) = u_complete(1:N);

%save data into simulation_params
simulation_params.u = u;
simulation_params.A = A;
simulation_params.B = B;
simulation_params.As = As;
simulation_params.s = s;
simulation_params.ns = ns;
simulation_params.cutoff = cutoff;