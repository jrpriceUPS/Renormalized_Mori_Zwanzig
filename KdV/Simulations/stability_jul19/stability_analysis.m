%a script to identify, for various values of epsilon, the smallest stable
%ROM

addpath ../../simulation_functions
addpath ../../nonlinear
addpath ../../analysis

clear all;close all;

epsilon = fliplr(0.065:0.005:0.1);
N_list = 8:2:36;

Re = (epsilon*2^(1/4)/(2*pi)).^-1;
stable_N = zeros(length(epsilon),1);

where_to_start = 0;

simulation_params.alpha = 1;      %coefficient on nonlinear term in KdV
simulation_params.dt = 1e-3;      %timestep
simulation_params.endtime = 100;   %end of simulation
simulation_params.howoften = 1;   %how often to save state vector
simulation_params.blowup = 1;     %if 1, instabilities cause simulation to end, but not give error
simulation_params.tol = 1e10;    %tolerance for identifying instabilities
simulation_params.initial_condition = @(x) sin(x);
simulation_params.name = 'complete';

correct_t_list_length = length(0:simulation_params.dt*simulation_params.howoften:simulation_params.endtime);

%run exact solution to time 10 and use to find t2 and t4 coefficients for
%each ROM
for j = 1:length(epsilon)
    min_stable = 0;
    i = where_to_start;
    simulation_params.epsilon = epsilon(j);
    
    while min_stable == 0
        i = i+1;
        N = N_list(i);
        [epsilon(j);N]
        simulation_params.N = N;          %number of positive modes to simulate
        
        [t_list,u_list] = KdV_solve(simulation_params);
        
        if length(t_list) == correct_t_list_length
           
            stable_N(j) = N
            min_stable = 1;
            where_to_start = i-1;
            
        elseif i == length(N_list)
            
            N_list(length(N_list)+1) = N_list(end)+2;
            
        end
    end
end

Lambda = 2*pi*stable_N;

Lambda_Re_ratio = Lambda./Re.'

save stable_N stable_N
save Lambda_Re_ratio Lambda_Re_ratio