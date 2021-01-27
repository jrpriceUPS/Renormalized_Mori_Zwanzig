function [simulation_params] = full_init_Burgers(simulation_params)
%
%[simulation_params] = full_init_Burgers(simulation_params)
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
%       N         =  number of positive resolved modes
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
%       M        =  size of the full model
%
%       u        =  state vector of positive modes of full model
%
%       RHS     =  function handle of the RHS

% Gather required parameters
N = simulation_params.N;
M = simulation_params.M;

% Define the ordinates in real space
x = linspace(0,2*pi*(2*N-1)/(2*N),2*N);

%define the initial condition as the Fourier transform of the sine function
%and ensure that high energy modes are zero in spite of any numerical
%errors
u_complete=fft_norm(simulation_params.initial_condition(x).');

%construct vector of modes we will advance in full model (just the positive modes), and
%fill it
u = zeros(N,1);
u(:) = u_complete(1:N);

%save data into simulation_params
simulation_params.u = u;
simulation_params.RHS = @(x,t) nonlinear_full_Burgers(x,t,simulation_params);