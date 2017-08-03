function simulation_params = complete_function_init(simulation_params)
%
%[simulation_params] = complete_fixed_init(simulation_params)
%
%Takes simulation_params structures and initializes them
%for a complete memory approximation with known effective coefficients simulation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
%
%       N         =  number of positive resolved modes
%
%       epsilon   =  degree of dispersion
%
%       alpha     =  coefficient of nonlinearity
%
%       dt        =  timestep
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
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
u_complete=fft_norm(sin(x).');

%initialize cells indicating index information, and populate them
simulation_params.F_modes = [1:N,2*N:4*N+2,5*N+2:6*N];
simulation_params.G_modes = N+1:5*N+1;
simulation_params.k = [0:3*N-1,-3*N:-1].';
simulation_params.M = 3*N;

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

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%construct vector of modes we will advance in full model (just the positive modes), and
%fill it
u = zeros(N,1);
u(:) = u_complete(1:N);

%save data into simulation_params
simulation_params.u = u;

simulation_params.coeffs = renormalized_ROM_coeffs(simulation_params.epsilon,simulation_params.N);

simulation_params.A = A;
simulation_params.B=@(x,t) renormalized_complete_4th_order_fixed_coeff(x,t,simulation_params);
simulation_params.As = As;
simulation_params.s = s;
simulation_params.ns = ns;
simulation_params.cutoff = cutoff;