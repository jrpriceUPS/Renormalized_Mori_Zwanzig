function simulation_params=nondim_BKdV(simulation_params)
%
%simulation_params=nondim_BKdV(simulation_params)
%
%A function to non-dimensionalize the units in a KdV-Burgers simulation

epsilon = simulation_params.epsilon;
U = sqrt(1/(2*pi)*integral(@(x) simulation_params.initial_condition(x).^2,0,2*pi));

x_scaling = epsilon/sqrt(U);
t_scaling = epsilon/(U^(3/2));

simulation_params.endtime = simulation_params.endtime / t_scaling;
simulation_params.dt = simulation_params.dt / t_scaling;
simulation_params.t_scaling = t_scaling;
simulation_params.alpha = epsilon*sqrt(U)/simulation_params.R;
simulation_params.L = 2*pi / x_scaling;
simulation_params.x_scaling = x_scaling;