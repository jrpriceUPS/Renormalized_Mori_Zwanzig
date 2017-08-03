function [renormalization_params,simulation_params] = renormalization_routine(renormalization_params,simulation_params)
%
%[renormalization_params,simulation_params] = renormalization_routine(renormalization_params,simulation_params)
%
%Takes initial renormalization_params and simulation_params structures and
%updates them according to the desired renormalization scheme.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
%       evolution         =  function handle for integration scheme for
%                            post-renormalization steps
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
%       F_modes  =  a cell array of indices for resolved modes in the full
%                   model
%
%       G_modes  =  a cell array of indices for unresolved modes in the full
%                   model
%
%       k        =  a cell array of wave numbers in full model
%
%       cutoff   =  cutoff for stiff portions of A
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% renormalization_params: a structure in which data relevant to the BCH
%                         renormalization scheme
%
%       RHS_match         =  updated vector of matching conditions from full model for
%                            each subset of modes used for renormalization
%
%       nonlinear_markov  =  updated vector of Markov contribution to energy flow
%                            for each subset of modes used for renormalization
%
%       nonlinear_terms   =  updated matrix of contribution to energy flow for each
%                            subset of modes used for renormalization due to
%                            each order of the BCH ROM expansion
%
%       done_yet          =  updated vector indicating whether data from each
%                            subset of modes has been saved for calculations
%
%       jettison          =  updated logical indicating whether the renormalization
%                            has been completed (1) or not (0)
%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
%
%       u       =  possibly updated state vector of positive modes of full model
%
%       A       =  possibly updated matrix of linear part of KdV model
%
%       B       =  possibly updated function handle of nonlinear part of RHS in full model
%
%       As      =  possibly updated stiff portion of A to be handled implicitly
%
%       s       =  possibly updated indices of stiff portion of A
%
%       ns      =  possibly updated indices of nonstiff portions of A
%
%       coeffs  =  computed renormalization coefficients

%create shorthands for ease of reading code
done_yet = renormalization_params.done_yet;
N_list = renormalization_params.N_list;
tol = renormalization_params.tol;

epsilon = simulation_params.epsilon;
alpha = simulation_params.alpha;
u = simulation_params.u;
cutoff = simulation_params.cutoff;

%loop through each mode subset we are matching
for j = 1:length(N_list)
    
    %if that mode has not been saved yet, inspect it
    if done_yet(j) ~= 1
        %identify mode of interest
        renormalization_params.index = j;
        N0 = N_list(j);
        
        %compute energy flow from subset for each ROM using the given
        %energy calculation scheme
        renormalization_params = renormalization_params.energy{renormalization_params.renorm_count+1}(simulation_params,renormalization_params);
        
        %if the t-model term passes the given threshold
        if norm(renormalization_params.nonlinear_terms(j,1)) > tol
            %mark this mode as completed
            renormalization_params.done_yet(j) = 1;
            %calculate the "full RHS" and save matching condition
            [RHS,~] = markov_term(u,3/2*renormalization_params.M,renormalization_params.M,alpha);
            renormalization_params.RHS_match(j) = sum(RHS(1:N0).*conj(u(1:N0))+conj(RHS(1:N0)).*u(1:N0));
            renormalization_params.RHS_match(j) = renormalization_params.RHS_match(j) - renormalization_params.markov_term(j);
        end
    end
end

%if the data for each subset of modes has been computed, calculate the
%renormalization coefficients
if sum(renormalization_params.done_yet) == length(N_list)
    %get rid of unresolved modes
    u = u(1:N_list(1));
    %calculate coefficients using the provided computation scheme
    coeffs = renormalization_params.coeff_scheme{renormalization_params.renorm_count+1}(renormalization_params.RHS_match,renormalization_params.nonlinear_terms,N_list)
    renormalization_params.renorm_count = renormalization_params.renorm_count+1;
    simulation_params.coeffs(:,renormalization_params.renorm_count) = coeffs;
    renormalization_params.done_yet(:) = 0;
    
    if renormalization_params.renorm_count == renormalization_params.total_renorm
        
        %mark the renormalization step as complete
        renormalization_params.jettison = 1;
        
        %update simulation_params in preparation for reduced model simulation
        simulation_params.N = N_list(1);
        simulation_params.F_modes = renormalization_params.F_modes_new;
        simulation_params.G_modes = renormalization_params.G_modes_new;
        simulation_params.k = renormalization_params.k_new;
        
        %update linear and nonlinear parts of RHS
        A = diag(1j*((1:N_list(1))-1).^3*epsilon^2);
        B = @(x,t) renormalization_params.evolution(x,t,simulation_params);
        
        if cutoff>N_list(1)
            ns  = 1:N_list(1);
            s = [];
        else
            ns  = 1:cutoff;
            s = cutoff+1:N_list(1);
        end
        As = A(s,s);
        
        %update simulation_params for reduced model simulation with new
        %evolution data
        simulation_params.u = u;
        simulation_params.A = A;
        simulation_params.B = B;
        simulation_params.As = As;
        simulation_params.s = s;
        simulation_params.ns = ns;
    end
end