function [t_list,u_list]=PDE_solve(simulation_params)
%
%[t_list,u_list]=PDE_solve(simulation_params,model)
%
%Solves a PDE using one of the many methods developed by myself
%and Panos.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% simulation_params: a structure in which data relevant to the simulation
%                    no matter the model type
%
%       epsilon    =  degree of dispersion
%
%       alpha      =  coefficient of nonlinearity
%
%       dt         =  timestep
%
%       starttime  =  time at which simulation begins (defaults to zero)
%
%       endtime    =  time at which simulation ends
%
%       howoften   =  how frequently to save the solution (i.e. 10 -> every 10
%                     timesteps)
%
%       blowup     =  logical variable indicating how to handle loss of energy
%                     conservation (if 1, save the data up to current time, if 0
%                     give error)
%
%       tol        =  tolerance for violation of conservation of energy
%
%       name       =  string indicating name of model (options: 'complete_fixed')
%
%       N          =  number of positive resolved modes
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Outputs:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  t_list  =  vector of the times at which data were recorded
%
%  u_list  =  an Nxlength(t_list) matrix containing the solution at each
%             time in t_list

%%%%%%%%%%%%%%%%
%Pre-processing%
%%%%%%%%%%%%%%%%

%general simulation parameters
dt       = simulation_params.dt;
endtime  = simulation_params.endtime;
howoften = simulation_params.howoften;
blowup   = simulation_params.blowup;
conservation_tolerance = simulation_params.tol;
N = simulation_params.N;

if ~isfield(simulation_params,'initial_condition')
    simulation_params.initial_condition = @(x) sin(x);
end

%specific model parameters
simulation_params = simulation_params.initialization(simulation_params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Time-Stepping Variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%construct vector of times
if ~isfield(simulation_params,'starttime')
    starttime = 0;
else
    starttime = simulation_params.starttime;
end
t = starttime:dt:endtime;
t_list = starttime:dt*howoften:endtime;

%identify the number of steps that need to be taken
nsteps = length(t)-1;

%save the initial total energy
total_energy = sum(2*abs(simulation_params.u).^2);

%create u storage matrix
u_list = zeros(N,length(t_list));
u_list(:,1) = simulation_params.u(1:N);
current_index = 1;

%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

for i=1:nsteps
    %update the current time
    simulation_params.current_time = t(i);
     
    %advance the simulation using the current ODE
    simulation_params.u = RK4_stiff_nonstiff_step(simulation_params);
    
    %if total energy in increases above the initial amount, abort the evaluation
    if ((sum(2*abs(simulation_params.u).^2)-total_energy)/total_energy>conservation_tolerance || isnan(sum(2*abs(simulation_params.u).^2)))
        
        if blowup
            u_list = u_list(:,1:current_index);
            t_list = t_list(1:current_index);
            break
        else
            error('Energy decay of t-model not preserved, try a smaller timestep')
        end
        
    end
    
    %save coefficients periodically
    if mod(i,howoften)==0
        current_index = current_index + 1;
        u_list(:,current_index) = simulation_params.u;
    end
end