function [t_list,u_list]=full_PDE_solve(simulation_params)
%
% [t_list,u_list]=full_PDE_solve(simulation_params,model)
%
% Solves for the full solution of PDE
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
%       starttime  =  time at which simulation begins (defaults to zero)
%
%       endtime    =  time at which simulation ends
%
%       name       =  string indicating name of model (options: 'complete_fixed')
%
%       N          =  number of positive resolved modes
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

options = odeset('AbsTol',1e-14,'Stats','on');
[t_list,u_list] = ode45(@(t,u) simulation_params.RHS(u,t),time_range,simulation_params.u,options);

u_list = permute(u_list,[2,1]);