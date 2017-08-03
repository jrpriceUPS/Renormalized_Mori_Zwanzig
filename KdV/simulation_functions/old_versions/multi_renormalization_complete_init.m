function [simulation_params,renormalization_params] = multi_renormalization_complete_init(model,simulation_params)
%
%[simulation_params,renormalization_params] = multi_renormalization_BCH_init(model,simulation_params)
%
%Takes initial model and simulation_params structures and initializes them
%for a complete renormalization simulation with two renormalization steps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% model: a structure indicating the model type and necessary
%        parameters
%
%  N             =  number of positive resolved modes
%
%
%  tol           =  tolerance for computing renormalization coefficients in
%                   renormalized models
%
%  coeff_scheme  =  mechanism for calculating renormalization constants
%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
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
%
%
% renormalization_params: a structure in which data relevant to the BCH
%                         renormalization scheme
%
%       N_list            =  vector of which subsets of modes to use for
%                            renormalization matching condition
%
%       RHS_match         =  vector of matching conditions from full model for
%                            each subset of modes used for renormalization
%
%       nonlinear_markov  =  vector of Markov contribution to energy flow
%                            for each subset of modes used for renormalization
%
%       nonlinear_terms   =  matrix of contribution to energy flow for each
%                            subset of modes used for renormalization due to
%                            each order of the BCH ROM expansion
%
%       done_yet          =  vector indicating whether data from each
%                            subset of modes has been saved for calculations
%
%       jettison          =  logical indicating whether the renormalization
%                            has been completed (1) or not (0)
%
%       tol               =  tolerance for norm of t-model energy drain at
%                            which that scales data is recorded
%
%       coeff_scheme      =  function handle for calculating renormalization 
%                            coefficients
%
%       energy            =  algorithm for calculating the energy flow from
%                            each subset of modes for each order of the BCH 
%                            ROM expansion

%create shorthand for ROM system size
N = simulation_params.N;

%general simulation parameters
epsilon = simulation_params.epsilon;
dt      = simulation_params.dt;
alpha   = simulation_params.alpha;


%define the ordinates in real space so there are 4N positive modes for full
%model
x=linspace(0,2*pi*(8*N-1)/(8*N),8*N);

%define the initial condition as the Fourier transform of the sine function
%and ensure that high energy modes are zero in spite of any numerical
%errors
u_complete=fft_norm(sin(x).');

%compute cutoff for different regimes of behavior
cutoff = ceil((2.8/(dt*epsilon^2))^(1/3));

if cutoff>4*N
    ns  = 1:4*N;
    s = [];
else
    ns  = 1:cutoff;
    s = cutoff+1:4*N;
end

%define the linear and nonlinear portions of the right hand side
A=diag(1j*((1:4*N)-1).^3*epsilon^2);
B=@(x,t) nonlinear_full(x,4*N,alpha);

%isolate the stiff portion of the linear portion of the right hand side
As = A(s,s);

%construct vector of modes we will advance in full model (just the positive modes), and
%fill it
u = zeros(4*N,1);
u(:) = u_complete(1:4*N);

%initialize cells indicating index information, and populate them
N_list = [N,N-1,N-2,N-3].';
F_modes = cell(length(N_list),1);
G_modes = cell(length(N_list),1);
k = cell(length(N_list),1);

for i = 1:length(N_list)
    %compute information about modes
    N0 = N_list(i);
    F_modes{i} = [1:N0,4*N0+1:4*N0+3,7*N0+4:8*N0+2];
    G_modes{i} = N0+1:7*N0+3;
    k{i} = [0:4*N0,-4*N0-1:-1].';
end


F_modes_new = [1:N,2*N+1:2*N+3,3*N+4:4*N+2];
G_modes_new = N+1:3*N+3;
k_new = [0:2*N,-2*N-1:-1].';

%save data into simulation_params
simulation_params.F_modes = F_modes;
simulation_params.G_modes = G_modes;
simulation_params.k = k;
simulation_params.u = u;
simulation_params.A = A;
simulation_params.B = B;
simulation_params.As = As;
simulation_params.s = s;
simulation_params.ns = ns;
simulation_params.cutoff = cutoff;
simulation_params.N = N_list(1);
simulation_params.M = 2*N_list(1)+1;

%save data into renormalization_params
renormalization_params.RHS_match = zeros(length(N_list),1);
renormalization_params.nonlinear_markov = zeros(length(N_list),1);
renormalization_params.nonlinear_terms = zeros(length(N_list),3);
renormalization_params.done_yet = zeros(length(N_list),1);
renormalization_params.jettison = 0;
renormalization_params.tol = model.tol;
renormalization_params.N_list = N_list;
renormalization_params.M = 4*N_list(1);
renormalization_params.coeff_scheme{1} = @(RHS_match,nonlinear_terms,N_list) coeffs_calculator(RHS_match,nonlinear_terms,N_list);
renormalization_params.coeff_scheme{2} = @(RHS_match,nonlinear_terms,N_list) coeffs_calculator_dispers(RHS_match,nonlinear_terms,N_list);
renormalization_params.energy{1} = @(a,b) multi_renormalized_complete_3rd_order_energy_nondispers(a,b);
renormalization_params.energy{2} = @(a,b) multi_renormalized_complete_3rd_order_energy_dispers(a,b);
renormalization_params.evolution = @(x,y,z) multi_renormalized_complete_3rd_order(x,y,z);
renormalization_params.renorm_count = 0;
renormalization_params.total_renorm = 2;
renormalization_params.F_modes_new = F_modes_new;
renormalization_params.G_modes_new = G_modes_new;
renormalization_params.k_new = k_new;