function [t_list,u_list]=ROM_PDE_solve(simulation_params)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialize Time-Stepping Variables%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = simulation_params.epsilon;

%construct vector of times
if ~isfield(simulation_params,'starttime')
    starttime = 0;
else
    starttime = simulation_params.starttime;
end


%%%%%%%%%%%%%%%%
%Pre-processing%
%%%%%%%%%%%%%%%%

if ~isfield(simulation_params,'initial_condition')
    simulation_params.initial_condition = @(x) epsilon*sin(x);
end

%specific model parameters
simulation_params = simulation_params.initialization(simulation_params);


if ~isfield(simulation_params,'time_range')
    time_range = [starttime,simulation_params.endtime];
else
    time_range = simulation_params.time_range;
end

%specific model parameters
simulation_params = simulation_params.initialization(simulation_params);


%%%%%%%%%%%%%%%%%%
%Integration Loop%
%%%%%%%%%%%%%%%%%%

percent_increase = 10;
init_energy = get_energy(simulation_params.u,simulation_params.N);
options = odeset('AbsTol',1e-14,'Stats','on','Events',@(t,u) energy_grow(t,u,init_energy,percent_increase,simulation_params));
[t_list,u_list] = ode45(@(t,u) simulation_params.RHS(u,t),time_range,simulation_params.u,options);

u_list = permute(u_list,[2,1]);